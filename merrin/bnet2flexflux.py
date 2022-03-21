
from xml.dom.minidom import parse, parseString

from colomoto.minibn import BooleanNetwork
import biolqm

note = parseString("""
    <notes>
        <body xmlns="http://www.w3.org/1999/xhtml">
            <p>STATE 0:[0,0]</p>
            <p>STATE 1:]0,+inf]</p>
        </body>
    </notes>
    """)
note = note.firstChild

def fixsbml(sbmlfile, ext, regulators):
    sbml = parse(sbmlfile)

    # compartments
    listNode = sbml.getElementsByTagName("listOfCompartments")[0]
    compDef = listNode.getElementsByTagName("compartment")[0]
    compDef.setAttribute("id", "cell")
    compExt = sbml.createElement("compartment")
    compExt.setAttribute("id", "ext")
    compExt.setAttribute("constant", "true")
    listNode.appendChild(compExt)
    for species in sbml.getElementsByTagName("qual:qualitativeSpecies"):
        name = species.getAttribute("qual:id")
        species.setAttribute("qual:compartment",
                "ext" if name in ext else "cell")
        species.setAttribute("qual:initialLevel", "0")
        if name not in regulators:
            species.appendChild(note.cloneNode(True))

    # remove transitions of ext species
    for qinput in sbml.getElementsByTagName("qual:output"):
        species = qinput.getAttribute("qual:qualitativeSpecies")
        if species in ext:
            qtr = qinput.parentNode.parentNode
            qtr.parentNode.removeChild(qtr)

    with open(sbmlfile, "w") as fp:
        sbml.writexml(fp)

def close_bn(bn):
    used = set()
    for f in bn.values():
        used.update([l.args[0].obj if isinstance(l, bn.ba.NOT) else l.obj
                for l in f.get_literals()])
    ext = used.difference(bn)
    for node in ext:
        bn[node] = node
    return ext

def bn2flexflux(bn, sbmlfile, regulators):
    close_bn(bn)
    print(bn)
    lm = bn.to_biolqm()
    biolqm.save(lm, sbmlfile)
    fixsbml(sbmlfile, bn.inputs() + list(bn.constants()), regulators)

def main():
    import sys
    bnetfile = sys.argv[1]
    sbmlfile = sys.argv[2]
    regulators = sys.argv[3:]
    bn = BooleanNetwork.load(bnetfile)
    bn2flexflux(bn, sbmlfile, regulators)

if __name__ == "__main__":
    main()
