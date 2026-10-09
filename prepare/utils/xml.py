"""The one lxml parser every PREPARE module reads an XML file with, and the function that uses it.

NCBI and store metadata files are mostly trusted, but they can be cached, replayed from a shared
store, or (for a store shared between projects) written by another process on the same host, so a
parser that would expand an entity or follow a DTD reference is a foothold for two well known XML
attacks: XXE (an external entity that reads a local file or calls out over the network) and
billion laughs (an internal entity that expands into gigabytes of text from a few bytes of
markup). ``SAFE_PARSER`` turns every one of those features off. The stdlib ``xml.etree`` call
sites elsewhere in PREPARE (``data/metadata.py``, ``data/sra_metadata.py``, ``data/taxonomy.py``)
are not moved onto this module: ElementTree never resolves an external entity (it raises an
undefined-entity ``ParseError``) or fetches a DTD, so XXE does not apply to those sites, and expat
2.4.1 and later rejects entity amplification such as billion laughs (its input amplification
limit). Those sites parse a small API response fetched directly from NCBI for immediate use, not a
file that could have been cached, replayed or written by another process.
"""

from pathlib import Path
from typing import Union

from lxml import etree

# resolve_entities=False is the main guard: an entity reference (external or internal,
# declared in the document's own internal DTD subset) is kept as an inert, unexpanded node
# instead of being substituted into the tree, which defeats XXE (the referenced file or URL is
# never read) and leaves a "billion laughs" document's entities unexpanded too, though libxml2's
# own entity-amplification limit (on regardless of these flags) is what actually raises on one
# nested deeply enough to matter. no_network=True refuses a DTD or entity fetched over the
# network even if resolve_entities were turned back on for one caller. load_dtd=False refuses an
# _external_ DTD or external parameter entity (a document's inline internal subset, the
# `<!DOCTYPE ... [ ... ]>` block, is still read to learn what entities exist; resolve_entities is
# what keeps their use inert). huge_tree=False (lxml's own default; named here for clarity rather
# than to change it) keeps libxml2's guards against a single grossly oversized document.
SAFE_PARSER = etree.XMLParser(resolve_entities=False, no_network=True, load_dtd=False, huge_tree=False)


def parse_xml_file(path: Union[str, Path]) -> etree._Element:
    """Parse the XML file at ``path`` with ``SAFE_PARSER`` and return its root element.

    Raises OSError when the file cannot be read and ``lxml.etree.XMLSyntaxError`` on malformed
    XML, the same errors a plain ``lxml.etree.parse(path).getroot()`` raises.
    """
    return etree.parse(str(path), parser=SAFE_PARSER).getroot()
