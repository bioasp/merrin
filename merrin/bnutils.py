def bn_score_influence_graph(bn, bn0):
    ig = bn.influence_graph()
    ig0 = bn0.influence_graph()
    e = set(((i,j,d["sign"]) for (i,j,d) in ig.edges(data=True)))
    e0 = set(((i,j,d["sign"]) for (i,j,d) in ig0.edges(data=True)))
    TP = e.intersection(e0)
    FP = e.difference(e0)
    return {
        "recall": len(TP)/len(e0),
        "precision": len(TP)/(len(TP)+len(FP)),
    }

if __name__ == "__main__":
    import sys
    from colomoto import minibn

    bn = minibn.BooleanNetwork.load(sys.argv[1])
    bn0 = minibn.BooleanNetwork.load(sys.argv[2])
    print(bn_score_influence_graph(bn, bn0))
