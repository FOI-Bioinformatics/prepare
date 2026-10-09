"""``prepare.utils.xml``: the shared lxml parser refuses entity expansion and network access.

Three attacks are checked directly against ``SAFE_PARSER``/``parse_xml_file`` rather than through
any one caller, since every lxml call site in PREPARE (currently ``data/metadata.py``) shares
this one parser: an internal entity expanded many times over ("billion laughs"), an external
entity reading a local file (XXE), and an external entity trying to reach the network.
"""

import time

import pytest
from lxml import etree

from prepare.utils.xml import SAFE_PARSER, parse_xml_file

# A small "billion laughs" document: 10 entities each expanding into 10 of the previous one, so
# the fully expanded text would be 10**9 repetitions of "lol" -- gigabytes from under 1 KB of
# markup, if the parser expanded it at all.
_BILLION_LAUGHS = """<?xml version="1.0"?>
<!DOCTYPE lolz [
 <!ENTITY lol0 "lol">
 <!ENTITY lol1 "&lol0;&lol0;&lol0;&lol0;&lol0;&lol0;&lol0;&lol0;&lol0;&lol0;">
 <!ENTITY lol2 "&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;&lol1;">
 <!ENTITY lol3 "&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;&lol2;">
 <!ENTITY lol4 "&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;&lol3;">
 <!ENTITY lol5 "&lol4;&lol4;&lol4;&lol4;&lol4;&lol4;&lol4;&lol4;&lol4;&lol4;">
 <!ENTITY lol6 "&lol5;&lol5;&lol5;&lol5;&lol5;&lol5;&lol5;&lol5;&lol5;&lol5;">
 <!ENTITY lol7 "&lol6;&lol6;&lol6;&lol6;&lol6;&lol6;&lol6;&lol6;&lol6;&lol6;">
 <!ENTITY lol8 "&lol7;&lol7;&lol7;&lol7;&lol7;&lol7;&lol7;&lol7;&lol7;&lol7;">
 <!ENTITY lol9 "&lol8;&lol8;&lol8;&lol8;&lol8;&lol8;&lol8;&lol8;&lol8;&lol8;">
]>
<root>&lol9;</root>
"""

_XXE_EXTERNAL_FILE = """<?xml version="1.0"?>
<!DOCTYPE root [
 <!ENTITY xxe SYSTEM "file:///etc/passwd">
]>
<root>&xxe;</root>
"""

_INTERNAL_ENTITY = """<?xml version="1.0"?>
<!DOCTYPE root [
 <!ENTITY greeting "hello">
]>
<root>&greeting;</root>
"""


def test_billion_laughs_is_rejected_quickly(tmp_path):
    """A billion-laughs document raises (libxml2's entity-amplification limit) in under 2 s, not a hang."""
    path = tmp_path / "billion_laughs.xml"
    path.write_text(_BILLION_LAUGHS, encoding="utf-8")

    started = time.monotonic()
    with pytest.raises(etree.XMLSyntaxError):
        parse_xml_file(path)
    elapsed = time.monotonic() - started

    assert elapsed < 2.0


def test_external_entity_is_not_expanded(tmp_path):
    """An external entity referencing a local file parses, but is never read or substituted.

    ``resolve_entities=False`` leaves the reference as an inert node instead of the file's
    content, so the joined text of the document is the literal ``&xxe;`` markup, never the
    target file's bytes; ``no_network=True`` is the same guard for a ``http://`` SYSTEM id.
    """
    path = tmp_path / "xxe.xml"
    path.write_text(_XXE_EXTERNAL_FILE, encoding="utf-8")

    root = parse_xml_file(path)

    text = "".join(root.itertext())
    assert text == "&xxe;"
    assert "root:" not in text  # a typical /etc/passwd line; never leaked into the tree


def test_internal_entity_is_left_unexpanded(tmp_path):
    """A harmless internal entity is left as literal markup too, not substituted into the tree.

    Contrasted with a parser that does resolve entities (``resolve_entities=True``), which
    substitutes ``&greeting;`` with ``hello`` as lxml's own default behaviour would.
    """
    path = tmp_path / "internal.xml"
    path.write_text(_INTERNAL_ENTITY, encoding="utf-8")

    root = parse_xml_file(path)
    text = "".join(root.itertext())
    assert text == "&greeting;"

    resolving_root = etree.fromstring(_INTERNAL_ENTITY.encode(), parser=etree.XMLParser(resolve_entities=True))
    assert resolving_root.text == "hello"


def test_well_formed_document_without_a_dtd_still_parses(tmp_path):
    """A normal NCBI-style document with no DOCTYPE parses exactly as before."""
    path = tmp_path / "plain.xml"
    path.write_text("<root><child>value</child></root>", encoding="utf-8")

    root = parse_xml_file(path)

    assert root.findtext("child") == "value"


def test_module_exposes_one_shared_parser_instance():
    """Every caller parses with the same ``SAFE_PARSER`` instance, not a fresh one per call."""
    assert isinstance(SAFE_PARSER, etree.XMLParser)


@pytest.mark.skipif(
    tuple(int(part) for part in __import__("pyexpat").version_info) < (2, 4, 1),
    reason="expat's amplification limit arrived in 2.4.1",
)
def test_stdlib_elementtree_rejects_billion_laughs_quickly():
    """The stdlib parse sites (NCBI responses) rely on expat's input amplification limit."""
    import xml.etree.ElementTree as ET

    started = time.monotonic()
    with pytest.raises(ET.ParseError):
        ET.fromstring(_BILLION_LAUGHS)
    assert time.monotonic() - started < 2.0


def test_stdlib_elementtree_does_not_resolve_an_external_entity(tmp_path):
    import xml.etree.ElementTree as ET

    secret = tmp_path / "secret.txt"
    secret.write_text("SECRET")
    document = f'<?xml version="1.0"?><!DOCTYPE r [<!ENTITY x SYSTEM "file://{secret}">]><r>&x;</r>'
    with pytest.raises(ET.ParseError, match="undefined entity"):
        ET.fromstring(document)
