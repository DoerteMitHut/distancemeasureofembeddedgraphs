
class Vertex():
    def __init__(self,index,x,y):
        self.index = index
        self.coords = (x,y)
        self.incidentEdges = set()
        self.placements = set()
        self.containingPlacements = dict()
        self.marked = False
        self.visited = False
        self.placementCount = 0;
        self.lowerBorder = 0;
        self.upperBorder = 1;
class Edge():
    def __init__(self,index,v1,v2):
        self.index = index
        self.v1 = v1
        self.v2 = v2
        self.marked = False
        self.containingPlacements = dict()
        self.edgePlacements = set()
        self.placementCount = 0
    def getOtherVertex(self,v):
        if self.v1 == v:
            return self.v2
        elif self.v2 == v:
            return self.v1
        else:
            raise Error("line " + str(self.index) + " has no Vertex with index " + str(v.getIndex())+"!")
    def getV1(self):
        return self.v1
    def getV2(self):
        return self.v2
    def getIndex(self):
        return self.index
    def getLength(self):
        return distance(self.v1,self.v2)
    def mark(self):
        self.marked = True
    def unmark(self):
        self.marked = False
    def isMarked(self):
        return self.marked
    def addContainingPlacement(self,placement,keyVertex):
        self.containingPlacements[keyVertex]=placement
    def addEdgePlacement(self,p):
        self.edgePlacements.add(p)
        p.setIndex(self.placementCount)
        self.placementCount = self.placementCount+1
    def getEdgePlacements(self):
        return self.edgePlacements
    def removeEdgePlacement(self,p):
        self.placementCount = self.placementCount-1
        self.edgePlacements.remove(p)
    def deleteAllEdgePlacements(self):
        self.marked = False
        self.containingPlacements = dict()
        self.edgePlacements = set()
        self.placementCount = 0
    def hasPlacements(self):
        if self.placementCount >0:
            return True
        else:
            return False
    def isEdgeInPlacementOf(self,v):
        if v in self.containingPlacements:
            return True
        else:
            return False

    def getContainingPlacement(self,keyVertex):
        return self.containingPlacements[keyVertex] 
class Placement():
    def __init__(self,u):
        self.index = 99999
        self.verts = set()
        self.edges = set()
        self.edgePlacements = set()
        self.valid = True
        self.marked = False
        self.degenerated = False
        self.motherVertex = u
    def getMotherVertex(self):
        return self.motherVertex
    def setDegenerated(self):
        self.degenerated = True
    def isDegenerated(self):
        return self.degenerated
    def getVerts(self):
        return self.verts
    def getEdges(self):
        return self.edges
    def getEdgePlacements(self):
        return self.edgePlacements
    def addVertex(self, v):
        self.verts.add(v)
    def addEdge(self,e):
        self.edges.add(e)
    def addEdgePlacement(self,p):
        self.edgePlacements.add(p)
    def removeEdgePlacement(self,ep):
        self.edgePlacements.remove(ep)
    def validate(self):
        self.valid = True
    def invalidate(self):
        self.valid = False
    def isValid(self):
        return self.valid
    def isMarked(self):
        return self.marked
    def mark(self):
        self.marked = True
    def unmark(self):
        self.marked = False
    def getIndex(self):
        return self.index
    def setIndex(self,index):
        self.index = index
    def __str__(self):
        width = 9+max(len(self.verts),len(self.edges))*3
        edgefill = width-(9+len(self.edges)*3)
        vertfill = width-(9+len(self.verts)*3)
        return (''.join(["Placement ",str(self.index)," of ",str(self.motherVertex.getIndex()),"\n","="*width,"\n|verts: ",str([v.getIndex() for v in self.verts])," "*vertfill,"|\n|edges: ",str([e.getIndex() for e in self.edges])," "*edgefill,"|\n","="*width]))
class EdgePlacement():
    def __init__(self,edge):
        self.index = 9999
        self.edge = edge
        self.v1Placements = set()
        self.v2Placements = set()
        self.v1Count = 0
        self.v2Count = 0
        self.valid = True
    def setIndex(self,i):
        self.index = i
    def getIndex(self):
        return self.index
    def addV1Placement(self,p):
        self.v1Placements.add(p)
        self.v1Count = len(self.v1Placements)
    def addV2Placement(self,p):
        self.v2Placements.add(p)
        self.v2Count = len(self.v2Placements)
    def getEdge(self):
        return self.edge
    def contains(self,p):
        if p in self.v1Placements or p in self.v2Placements:
            return True
        else:
            return False
    def getOtherSide(self,p):
        if p in self.v1Placements:
            return list(self.v2Placements)
        elif p in self.v2Placements:
            return list(self.v1Placements)
        else:
            raise Error("edgePlacement of edge" +str(self.edge().getIndex())+ " does not contain the placement you want to remove " + str(v.getIndex())+"!")
    def getV1(self):
        return list(self.v1Placements)
    def getV2(self):
        return list(self.v2Placements)
    def validate(self):
        self.valid = True
    def invalidate(self):
        self.valid = False
    def isValid(self):
        return self.valid
    def getV1Count(self):
        return self.v1Count
    def getV2Count(self):
        return self.v2Count
    def removePlacement(self,p):
        if p in self.v1Placements:
            self.v1Placements.remove(p)
            self.v1Count = self.v1Count-1
        elif p in self.v2Placements:
            self.v2Placements.remove(p)
            self.v2Count = self.v2Count-1 
        else:
            raise Error("edgePlacement of edge" +str(self.edge.getIndex())+ " does not contain the placement you want to remove " + str(v.getIndex())+"!")
        if self.v1Count <= 0 or self.v2Count <= 0:
            self.invalidate()
    def __str__(self):
        return(''.join(["placement of edge ",str(self.edge.getIndex()),":\n\tPlacements of vertex ",str(self.edge.getV1().getIndex())," are:\n\t",str([p.getIndex() for p in self.v1Placements]),"\n\tPlacements of vertex ",str(self.edge.getV2().getIndex())," are:\n\t",str([p.getIndex() for p in self.v2Placements])]))
class graphRep():
    def __init__(self):
        self.verts = {}
        self.edges = {} 
        self.coordinates = []*0
        self.hasVerts = False
        self.hasEdges = False
        self.vertFile = ""
        self.edgeFile = ""

def mappingBordersExp(v,e,epsilon,p):
    #computes the explicit borders for a mapping of vertex v onto edge e within a ball of radius epsilon.
    startvert = p.getMotherVertex()
    #incident vertex to e, that marks position 0 on the edge.
    destvert = e.getOtherVertex(startvert)
    #incident vertex to e, that marks position 1 on the edge.
    anchor = perpendicularFoot(v,e)
    #compute the orthogonal projection of v onto e
    d = distance(v,anchor)
    #compute length of projection vector
    if d > epsilon:
        #the verte is too far away from the edge and cannot be mapped.
        return[1,0]#CHECK
    offset = sqrt(pow(max(d,epsilon),2)-pow(min(d,epsilon),2))
    # max(d,epsilon) describes the hypothenuse of a Triangle between the embedding of v, its projection onto e and
    # the intersection of the epsilon circle around v with the embedding of e.
    # min(d,epsilon) describes the cathede which does not lie inside of e. Offset therefore values at half the length of the final mapping.
    nVec = [(destvert.getCoords()[0]-startvert.getCoords()[0])/e.getLength(),(destvert.getCoords()[1]-startvert.getCoords()[1])/e.getLength()]
    # normalized vector parallel to e 
    hiMap = Vertex(99999,anchor.getCoords()[0]+offset*nVec[0],anchor.getCoords()[1]+offset*nVec[1])
    # construct vertex at higher mapping limit
    loMap = Vertex(99999,anchor.getCoords()[0]-offset*nVec[0],anchor.getCoords()[1]-offset*nVec[1])
    # construct vertex at lower mapping limit
    if distance(startvert,hiMap)>=distance(destvert,hiMap):
        # destvert closer to hiMap or hiMap at center of e
        if distance(startvert,hiMap)>= e.getLength():
            # hiMap outside of e beyond destvert
            hi = 1
        else:
            # hiMap in the destvert half of e
            hi = 1-distance(destvert,hiMap)/e.getLength()
    else:
        # hiMap closer to startvert
        if distance(destvert,hiMap)>e.getLength():
            # mapping outside of e
            return [1,0]#CHECK //hi = 0-distance(startvert,hiMap)/e.getLength()
        else:
            # hiMap in the startvert haf of e
            hi = 0+distance(startvert,hiMap)/e.getLength()

    if distance(destvert,loMap)>=distance(startvert,loMap):
        # loMap closer to startvert
        if distance(destvert,loMap)>= e.getLength():
            # loMap outside of e beyond startvert
            lo = 0
        else:
            # loMap inside startvert half of e
            lo = 0+distance(startvert,loMap)/e.getLength()
    else:
        # loMap closer to destvert
        if distance(startvert,loMap)>e.getLength():
            # loMap outside e beyond destvert. mapping outside e
            return [1,0]#CHECK //lo = 1+distance(destvert,loMap)/e.getLength()
        else:
            # loMap in the destvert half of e
            lo = 1-distance(destvert,loMap)/e.getLength()
    return [lo,hi] 





def findPlacement(u, v, epsilon):
    # begins to traverse the target graph beginning from v.
    # constructs placement of u in g by traversing unmarked vertices of g
    # inside the epsilon-ball
    if distance(u,v) <= epsilon:    
    #first vertex is checked for being inside of the epsilon-ball
        stack = []                  
        # a stack is initialized to serve the traversal
        p = Placement(u)
        # a placement of u in g is initialized
        p.addVertex(v)
        # v is added to the placement
        u.addPlacement(p)
        # the placement and its mother vertex are linked
        v.addContainingPlacement(p,u)
        # placement components store their containing placements 
        #keyed with the mother vertices
        v.mark()
        # the first vertex of placement p is marked as visited

        for e in v.getEdges():      
            # all incident edges are added to the placement
            p.addEdge(e)
            e.mark()
            # the added edges are marked
            stack.append(e.getOtherVertex(v))
            # all adjacent vertices are pushed to the stack
        while len(stack)>=1:
            # the stack is processed
            currentVertex = stack.pop()
            if not currentVertex.isMarked():
                # marked vertices occur as parts of already existing
                # placements or inside of graphs with cycles
                # they are skipped
                if distance(u,currentVertex) <= epsilon:
                    # unmarked vertices within the epsilon-ball...
                    currentVertex.addContainingPlacement(p,u)
                    p.addVertex(currentVertex)
                    # ...are added to placement p
                    for e in currentVertex.getEdges():
                        # all incident edges are added to the placement,
                        p.addEdge(e)
                        e.addContainingPlacement(p,u)
                        e.mark()
                        # marked,
                        stack.append(e.getOtherVertex(currentVertex))
                        # and their 2nd vertex is pushed down the stack
                currentVertex.mark()
    else:
        # if v is not close enough to u,
        # it is marked for not being part of any placement of u.
        v.mark()
def findEdgePlacement(p,c,epsilon):
    # constructs a placement for an edge of f from 
    # a scaffolding c with one starting vertex placement p
    stack = list(random.sample(p.getVerts(),1))
    # initialize stack with a random vertex from the starting placement
    if p.isDegenerated() and not inEpsilonTube(stack[0],c.getEdge(),epsilon):
        # in order for the starting vertex placement of one incident vertex of the to-be-placed edge
        # to be a valid part of the edge placement in construction,
        # it has to be connected to a placement of the edges' other vertex by a path inside the epsilon-tube.
        # As regular vertex placements contain only vertices inside the epsilon ball around the mother vertex,
        # it doesn't matter with which one to begin traversal.
        # As far as degenerated placements go, they contain only one edge, the incident vertices of which are part of
        # the vertex placement, without actually being situated within the epsilon-ball. Should both incident vertices lie
        # outside the tube, the vertex placement would be an invalid part of the edge placement in construction.
        # Should one of the vertices lie outside, the vertex placement might illegally get folded, if we randomly begin at the wrong vertex.
        # we therefore have to try both vertices to handle this special case. 
        v = stack.pop()
        stack.append(list(p.getEdges())[0].getOtherVertex(v))
        # replace wrong starting vertex with correct one
    while len(stack) >= 1:
        # begin traversal inside of epsilon-tube
        currentVertex = stack.pop()
        if not currentVertex.isMarked():
            # only unvisited vertices are handled
            currentVertex.mark()
            vertexGood = False
            # each vertex has to falsify its' invalidity by being inside the epsilon-tube 
            if currentVertex.isInPlacementOf(c.getEdge().getV1()):
                # vertices encountered by the traversal might be part of a placement of the edges' first vertex.
                v1p = currentVertex.getContainingPlacement(c.getEdge().getV1())
                # as each vertex of g may only be part of one placement for each of the vertices in f,
                # each one of a g-vertices' containing placements can be sufficiently identified by it's mother vertex from f. 
                # v1p is a potential member of the v1-side of the edge placement in construction.
                if inEpsilonTube(currentVertex,c.getEdge(),epsilon):
                    # the current vertex has to be inside the epsilon-tube around the mother edge of the edge placement in construction.
                    c.addV1Placement(v1p)
                    # the placement is connected to the starting placement withon the epsilon-tube and joins the v1-side of the edge placement in construction.
                    vertexGood = True
                else:# I think this special cas is not special and does not need custom handling
                    # # if the encountered vertex is part of a placement of v1 but lies outside the epsilon tube...
                    # if v1p.isDegenerated():
                    #     # and if the placement in question is also degenerated, the other vertex might still be part of a placement. 
                    #     if inEpsilonTube(list(v1p.getEdges())[0].getOtherVertex(currentVertex),c.getEdge(),epsilon):
                    #         c.addV1Placement(v1p)
                    #         currentVertex = list(v1p.getEdges())[0].getOtherVertex(currentVertex)
                    #         currentVertex.mark()
                    #         vertexGood = True
            if currentVertex.isInPlacementOf(c.getEdge().getV2()):
                # The current vertex might also be part of a placement of the second vertex of the mother edge, V2.
                v2p = currentVertex.getContainingPlacement(c.getEdge().getV2())
                # the containing placement is appropriately queried
                if inEpsilonTube(currentVertex,c.getEdge(),epsilon):
                    # if the vertex is actually reachable within the epsilon-tube, it joins the edge placements' v2-side
                    c.addV2Placement(v2p)
                    vertexGood = True
                else:# strange special case like in the v1-section
                    # if v2p.isDegenerated():
                    #     if inEpsilonTube(list(v2p.getEdges())[0].getOtherVertex(currentVertex),c.getEdge(),epsilon):
                    #         c.addV1Placement(v2p)
                    #         currentVertex = list(v2p.getEdges())[0].getOtherVertex(currentVertex)
                    #         currentVertex.mark()
                    #         vertexGood = True
            else:
                # even if the vertex is not part of either placement, it might be situated inside the epsilon tube 
                # and marked as part of the mapped path
                if inEpsilonTube(currentVertex,c.getEdge(),epsilon):
                    vertexGood = True
                    currentVertex.mark()
            if vertexGood:
                #adjacent vertices are only added, if the current vertex is inside the epsilon tube.
                stack = stack+[e.getOtherVertex(currentVertex) for e in currentVertex.getEdges()]
def findEdgePlacementStrong(p,c,epsilon):
    lastVert = random.sample(p.getVerts(),1)[0]
    # get random vertex from starting vertex placement
    if p.isDegenerated() and not inEpsilonTube(lastVert,c.getEdge(),epsilon):
        # switch vertices if placement is degenerated and initial vertex is not inside the epsilon tube. 
         lastVert = list(p.getEdges())[0].getOtherVertex(lastVert)
         if not inEpsilonTube(lastVert,c.getEdge(),epsilon):
            # discard edge placement if both vertices of a degenerated starting placement are not inside the epsilon tube
            #CHECK: added double degenerated one-edge edge placements
            degE = list(p.getEdges())[0]
            if degE.isEdgeInPlacementOf(degE.getOtherVertex(lastVert)) and degE.getContainingPlacement(degE.getOtherVertex(lastVert)).isDegenerated():
                if p.getMotherVertex() == c.getEdge().getV1():
                    c.addV2Placement(degE.getContainingPlacement(degE.getOtherVertex(lastVert)))
                    c.addV1Placement(p.getMotherVertex())
                else:
                    c.addV1Placement(degE.getContainingPlacement(degE.getOtherVertex(lastVert)))
                    c.addV2Placement(p.getMotherVertex())
                degE.getContainingPlacement(degE.getOtherVertex(lastVert)).addEdgePlacement(c)
                p.getMotherVertex().addEdgePlacement(c)


    lastVert.setBorders(mappingBordersExp(lastVert,c.getEdge(),epsilon,p))
    # set borders by calling method to calculate explicit borders.
    stack = []
    #initialize empty stack
    if lastVert.isInPlacementOf(c.getEdge().getOtherVertex(p.getMotherVertex())):
        # check, wether the last visited vertex is in a placement of the destination vertex of e
        v2p = lastVert.getContainingPlacement(c.getEdge().getOtherVertex(p.getMotherVertex()))
        # get this placement
        #v2p.mark()
        if c.getEdge().getOtherVertex(p.getMotherVertex()) == c.getEdge().getV1():
            # if the last vertex is part of a placement of the vertex incident to c's mother edge, that the traversal wasn't started from,
            # it is either e's V1 or V2. It is added to the according side of the edge placement.
            c.addV1Placement(v2p)
            v2p.addEdgePlacement(c)
            #add V1 Placements to c tat contain currentVertex
        elif c.getEdge().getOtherVertex(p.getMotherVertex()) == c.getEdge().getV2():
            # other case. lastvertex is part of V2's placement.
            c.addV2Placement(v2p)
            v2p.addEdgePlacement(c)
        else:
            # as "isInPlacementOf" was True, this cannot be and indicates an ERROR.
            raise Error("FEHLER")

    for v in [e.getOtherVertex(lastVert) for e in lastVert.getEdges()]:
    # iterate over all adjacent vertices of lastVert
    #CHECK: this comment looks scary but I think the computation of borders returns an invalid mapping for everything outside the epsilon-tube.
    #if inEpsilonTube(v,c.getEdge(),epsilon):
        v.setBorders(mappingBordersExp(v,c.getEdge(),epsilon,p))
        # v gets it's explicit borders set.
        v.setLowerBorder(max(v.getLowerBorder(),lastVert.getLowerBorder()))
        # its actual lower border is the maximum of the explicit one and the preceding vertex's lower border.
        #CHECK: As this always maps preceding vertices first and restricts succeeding vertices, this might compromise the validity of computing just one direction. 
        if v.getUpperBorder()>=v.getLowerBorder():
            # there is a mapping for this vertex.
            stack.append(v)
        else:
            # there is no mapping this vertex onto the edge and it it branded accordingly.
            v.mark()
    #else:
        #print("...not in tube")
        #v.mark()
    lastVert.mark()
    # lastVert and it's adjacent vertices have been updated and stacked as far as their borders are concerned.
    while len(stack) >= 1:
        # clear the stack
        currentVertex = stack.pop()
        if not currentVertex.isMarked():
            #consider only unmarked vertices
            if currentVertex.isInPlacementOf(c.getEdge().getOtherVertex(p.getMotherVertex())):
                v2p = currentVertex.getContainingPlacement(c.getEdge().getOtherVertex(p.getMotherVertex()))
                if c.getEdge().getOtherVertex(p.getMotherVertex()) == c.getEdge().getV1():
                    c.addV1Placement(v2p)#add V1 Placements to c tat contain currentVertex
                    v2p.addEdgePlacement(c)
                elif c.getEdge().getOtherVertex(p.getMotherVertex()) == c.getEdge().getV2():
                    c.addV2Placement(v2p)
                    v2p.addEdgePlacement(c)
                else:
                    raise Error("FEHLER")
            for v in [e.getOtherVertex(currentVertex) for e in currentVertex.getEdges()]:
                v.setBorders(mappingBordersExp(v,c.getEdge(),epsilon,p))
                v.setLowerBorder(max(v.getLowerBorder(),currentVertex.getLowerBorder()))
                if v.getUpperBorder()>=v.getLowerBorder():
                    stack = stack + [v]
            currentVertex.mark()
def constructVertexPlacements(f,g,epsilon):
    # constructs vertex placements for the graph representations f and g
    # with proposed maximum distance epsilon 
    for u in f.getVerts():
        #iterate over all vertices of f
        for v in g.getVerts():
            #iterate over all vertices of g
            if not v.isMarked():
                # vertices from g may already be marked.
                # if this is not the case, a new placement
                # will be constructed via traversal 
                findPlacement(u,v,epsilon)
                # an ordinary placement is constructed or filed,
                # if v is not inside the epsilon-ball around v.
        for e in g.getEdges():
            # degenerated placements call for iterating over all edges of g
            if not e.isMarked():
                # if e is not already part of a placement of u,
                findDegPlacement(u,e,epsilon)
                # it is checked for being a degenerated placement
        if not u.hasPlacements():
            # abort if there are no epsilon-placements for u in g
            return False
        g.unmarkEverything()
        # else unmark everything in g to contruct placements
        # for the next vertex of f 

def constructEdgePlacements(f,g,epsilon):  
    for e in f.getEdges():
        #iterate over all edges of f
        for p in e.getV1().getPlacements():
            # iterate over placements of 1st vertex of current edge
            if not p.isMarked():
                # only unmarked vertex placements are considered for
                # beginning construction of the current edge placement
                c = EdgePlacement(e)
                # initialize new edge placement
                c.addV1Placement(p)
                # add current vertex placement to edge placement in construction
                findEdgePlacement(p,c,epsilon)
                # construct edge placement
                if c.getV1Count() != 0 and c.getV2Count() != 0:
                    # an edge placement has to have at least one placement of either incident vertex of its' mother edge in order to be valid.
                    e.addEdgePlacement(c)
                    # the mother edge and the constructed edge placement are linked
                    for vp in c.getV1():
                        # the constructed edge placement and its' v1 placements are linked 
                        vp.addEdgePlacement(c)
                    for vp in c.getV2():
                        # the constructed edge placement and its' v1 placements are linked
                        vp.addEdgePlacement(c)
                p.mark()
        f.unmarkEverything()
        g.unmarkEverything()
def constructEdgePlacementsStrong(f,g,epsilon):
    for e in f.getEdges():
        #iterate over all edges of f
        for p in e.getV1().getPlacements():
            # get all v1 placements
            c = EdgePlacement(e)
            # initialize edge placement
            c.addV1Placement(p)
            p.addEdgePlacement(c)
            # link with one vertex placement of v1
            findEdgePlacementStrong(p,c,epsilon)
            # construct edge placement c
            if c.getV1Count() != 0 and c.getV2Count() != 0:
                # neither side of edge placement c is empty
                e.addEdgePlacement(c)
                # add placement
            else:
                # edge placement has at least one empty side
                p.removeEdgePlacement(c)
                # unlink p and c
                #CHECK: there there were probably placements added to the non-empty side, which are definetely invalid, but are still linked to c.
                del c
            f.unmarkEverything()
            g.unmarkEverything()
        for p in e.getV2().getPlacements():
            #CHECK: as all V1 placements which are connected inside the tube were either used as starting points for a traversal or have been discovered during one, 
            # starting traversals from the V2 side should be obsolete.
            c = EdgePlacement(e)
            c.addV2Placement(p)
            p.addEdgePlacement(c)
            findEdgePlacementStrong(p,c,epsilon)
            if c.getV1Count() != 0 and c.getV2Count() != 0:
                e.addEdgePlacement(c)
            else:
                p.removeEdgePlacement(c)
                del c
        #==========================================================================================
        # for degE in g.getEdges():
        #     # traverse all edges of g
        #     if degE.isEdgeInPlacementOf(e.getV1()) and degE.isEdgeInPlacementOf(e.getV2()):
        #         # having a single edge cross both epsilon balls constructs a degenerated edge placement which hasn't been visited yet.
        #         # This is because neither one of the incident vertices has a mapping on e. 
        #         # CHECK: This spawns new edge placements all edges, which intersect with both epsilon balls. Those with incident vertices
        #         # inside the epsilon balls are also considered, despite them having already been counted towards the regular placements. 
        #         c = EdgePlacement(e)
        #         c.addV1Placement(degE.getContainingPlacement(e.getV1()))
        #         degE.getContainingPlacement(e.getV1()).addEdgePlacement(c)
        #         c.addV2Placement(degE.getContainingPlacement(e.getV2()))
        #         degE.getContainingPlacement(e.getV2()).addEdgePlacement(c)
        #         e.addEdgePlacement(c)
        #==========================================================================================
            f.unmarkEverything()
            g.unmarkEverything()

def pruneWeak(epsilon,f):
    # In order to determine a manifold of mappings between f and g, we prune away placements
    # which do not have paths to at least one placement of each of their mother vertex's adjecent vertices.
    stack =set([p for v in f.getVerts() for p in v.getPlacements()])
    # initialize a set of all vertex placements.
    while len(stack)>=1:
        # one placement at a time
        currentPlacement = stack.pop()
        # evaluate validity
        if currentPlacement.isValid():
            # for placements that have not been invalidated, check all edges incident to their mother vertices 
            for e in currentPlacement.getMotherVertex().getEdges():
                # edges are invalid by default
                eValid = False
                # check all placements of the current edge
                # at least one placement of each edge has to contain the current placement.
                #==========================================================================
                # maybe query the containing placements of the mother vertex with O(n^2) possible incident edges instead.
                # currently each vertex has to be checked for O(n^2) edges with O(n^2) possible edge placements each.
                # Seems not so smart.
                #==========================================================================
                for ep in e.getEdgePlacements():
                    if ep.contains(currentPlacement):
                        eValid = True
                        break
                if not eValid:
                    # if an edge turns out invalid, the vertex placement is as well and it is invalidated.
                    currentPlacement.invalidate()
                    break
            if not currentPlacement.isValid():
                # if the placement has been invalidated, its' invalidity propagates to all edge placements which rely on it as their only representative of one side.
                # invalidated edge placements discontinue their support for incident vertex placements, which in turn may invalidate those.
                currentPlacement.getMotherVertex().removePlacement(currentPlacement)
                # the link between the placement and its' mother vertex is destroyed
                for ep in currentPlacement.getEdgePlacements().copy():
                    # iterate over a list of all edge placements containing the invalidated vertex placement
                    ep.removePlacement(currentPlacement)
                    currentPlacement.removeEdgePlacement(ep)
                    # destroy their mutual link
                    if not ep.isValid():
                        # the edge placement may have been invalidated in this process. if this is the case:
                        ep.getEdge().removeEdgePlacement(ep)
                        # destroy the link between the mother edge and the edge placement
                        for p in ep.getV1():
                            # remove all v1 placements from the edge placement 
                            p.removeEdgePlacement(ep)
                            ep.removePlacement(p)
                            # As each vertex placement appears at most in one edge placement per edge incident to its' mother vertex,
                            # the destruction of an edge placement invalidates all contained vertex placements [REVIEW]
                            p.invalidate()
                        for p in ep.getV2():
                            # same as above, but for v2 placements
                            p.removeEdgePlacement(ep)
                            ep.removePlacement(p)
                            p.invalidate()
                        del ep # the emptied edge placement is discarded
    toRet = True
    # The pruning reports the decision algorithm's outcome. It is "True", unless a vertex of f lost all its' vertex placements.
    for v in f.getVerts():
        if not v.hasPlacements():
            toRet = False
            break
    return toRet
def pruneStrong(epsilon,f):
    stack =set([p for v in f.getVerts() for p in v.getPlacements()])
    while len(stack)>=1:
        currentPlacement = stack.pop()
        if currentPlacement.isValid():
            for e in currentPlacement.getMotherVertex().getEdges():
                eValid = False
                for ep in e.getEdgePlacements():
                    if ep.contains(currentPlacement):
                        eValid = True
                        break
                if not eValid:
                    currentPlacement.invalidate()
                    break
            if not currentPlacement.isValid():
                currentPlacement.getMotherVertex().removePlacement(currentPlacement)
                for ep in currentPlacement.getEdgePlacements().copy():
                    ep.removePlacement(currentPlacement)
                    currentPlacement.removeEdgePlacement(ep)
                    if not ep.isValid():
                        ep.getEdge().removeEdgePlacement(ep)
                        for p in ep.getV1():
                            p.removeEdgePlacement(ep)
                            stack.add(p)
                            ep.removePlacement(p)
                        for p in ep.getV2():
                            p.removeEdgePlacement(ep)
                            stack.add(p)
                            ep.removePlacement(p)
                        del ep
    toRet = True
    for v in f.getVerts():
        if not v.hasPlacements():
            toRet = False
            break
    return toRet
def decision(f,g,epsilon,strong):
    # f and g are Graph representations
    #
    # epsilon is the proposed distance, that has to be subceeded 
    # in order for the function to return "True" 
    #
    # strong is a boolean that states, whether the ordinary or 
    # "weak" frechet distance is to be used     
    vpPossible = constructVertexPlacements(f,g,epsilon)
    if vpPossible:
        if not strong:
            constructEdgePlacements(f,g,epsilon)
            return pruneWeak(epsilon,f)
        else:
            constructEdgePlacementsStrong(f,g,epsilon)
            return pruneStrong(epsilon,f)
    else:
        return Fals