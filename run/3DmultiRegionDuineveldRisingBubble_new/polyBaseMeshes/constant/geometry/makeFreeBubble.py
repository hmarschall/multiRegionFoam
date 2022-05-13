# This script creates a mesh of cubic domain. In the center a bubble
# is located, sourrounded by a fluid, both parts are meshed with
# tetrahedrons. Both boundary shells at the interface are meshed with
# radial prisms.

def makeBubble():
	import salome
	salome.salome_init()

	import GEOM
	
	from salome.geom import geomBuilder
	geompy = geomBuilder.New(salome.myStudy)

	import SMESH
	from salome.smesh import smeshBuilder
	smesh = smeshBuilder.New(salome.myStudy)
	
	import math
	from math import sqrt,pow

	# -------------------------------------------------------------------- #
	# Settings
	# -------------------------------------------------------------------- #


	# Radius of the bubble (in mm)
	# All other quantities are relative to the bubble radius
	bubbleRadius			= 0.001 #section

	# Ratio of the radius of the whole domain and the bubble radius
	relFluidRadius			= 20

	# Path and filename where the mesh shall be saved as UNV
	savePathUNV = '/home/local/CSI/tp83ytyb/foam/tp83ytyb-3.1-sfb1194/run/reactiveBubble/smallTestBubble/base/constant/geometry/freeBubble.unv'

	# Path and filename where the mesh shall be saved as STL
	savePathSTL = '/home/local/CSI/tp83ytyb/foam/tp83ytyb-3.1-sfb1194/run/reactiveBubble/smallTestBubble/base/constant/geometry/freeBubble.stl'

	# Approximate number of faces that bubble and liquid should have
	numFacesBubble = 20000 #12000 
	numFacesLiquidIO = 100000 #500 
	#numFacesLiquidW = 1250 #1250	

	# Approximate number of cells of bubble and liquid
	numElemsBubble = 10000 #12000 
	numElemsLiquid = 100000 #80000
	

	# If equal to 1 the surface meshing will produce also quadrangles,
	# but exporting to UNV is then impossible
	allowQuadrangles = 0
	
	# -------------------------------------------------------------------- #
	# Debug params
	# -------------------------------------------------------------------- #

	# Ratio of gap between layers 2 and 3 (non-debug=0)
	debug_gap	= 0



	# -------------------------------------------------------------------- #
	# Geometry
	# -------------------------------------------------------------------- #
	
	# Length of box edge (in mm)
	#rBox = 1.2*bubbleRadius 
	rSphere = 20*bubbleRadius 
	
	#Create sphere as revolution of semicircle

	Vertex_1 = geompy.MakeVertex(0, 0, 0)
	Vertex_2 = geompy.MakeVertex(0, 0, 0.01)
	Vertex_3 = geompy.MakeVertex(0, bubbleRadius, 0)
	Vertex_4 = geompy.MakeVertex(0, 0, -bubbleRadius)
	Vertex_5 = geompy.MakeVertex(0, 0, bubbleRadius)
	vector = geompy.MakeVector(Vertex_1, Vertex_2)
	
	points = []
	points.append(Vertex_1)
	points.append(Vertex_2)
	points.append(Vertex_3)

	#smoothingsurface = geompy.MakeSmoothingSurface( points )
	#id_smoothingsurface = geompy.addToStudy(smoothingsurface,"SmoothingSurface")

	circle = geompy.MakeArc(Vertex_4, Vertex_3, Vertex_5)
	geompy.addToStudy(circle, "Circle")
	
	#vector = geompy.MakeVectorDXDYDZ(0, 0, 0.001)
	geompy.addToStudy(vector, "revolutionAxis")
	freeBubble = geompy.MakeRevolution(circle, vector, 360*math.pi/180.0)
	geompy.addToStudy(freeBubble, "freeBubble")

	
	# Create external cylinder
	sphere = geompy.MakeSphereR(rSphere)
	geompy.addToStudy(sphere, "sphere")

	# Create the shells for the Taylor Bubble and the channel

	freeBubbleShell = geompy.MakeShell([freeBubble])
	freeBubbleSolidShadow = geompy.MakeSolid([freeBubbleShell])
	freeBubbleSolid = geompy.MakeSolid([freeBubbleShell])

	shells = []
	shells.append(freeBubbleSolidShadow)
	shells.append(geompy.MakeCut(sphere,freeBubbleSolid))

	# Create partitions and compound
	innerPartition = geompy.MakePartition([shells[0]])
	outerPartition = geompy.MakePartition([shells[1]])
	domainCompound = geompy.MakeCompound([innerPartition,outerPartition])
	id_domainCompound = geompy.addToStudy(domainCompound,"domainCompound")

	# Get faces
	# Faces 2 and 3 need to be swapped since salome changed order
	faces=geompy.SubShapeAll(domainCompound,geompy.ShapeType["FACE"])
	swap=faces[2]
	faces[2]=faces[1]
	faces[1]=swap
	for faceI,face in enumerate(faces):
		name = "Face {0}".format(faceI)
		geompy.addToStudyInFather(domainCompound, face, name)
		
	shells=geompy.SubShapeAll(domainCompound,geompy.ShapeType["SOLID"])
	for shellI,shell in enumerate(shells):
		name = "Shell {0}".format(shellI)
		geompy.addToStudyInFather(domainCompound, shell, name)

	
	# Create the groups
	groups = []
	
	# interfaceShadow group
	groups.append(geompy.CreateGroup(domainCompound, geompy.ShapeType["FACE"]))
	FaceID = geompy.GetSubShapeID(domainCompound, faces[0])
	geompy.AddObject(groups[0], FaceID)
	geompy.addToStudyInFather(domainCompound,groups[0], "Group interfaceShadow")
	
	# interface group
	groups.append(geompy.CreateGroup(domainCompound, geompy.ShapeType["FACE"]))
	FaceID = geompy.GetSubShapeID(domainCompound, faces[1])
	geompy.AddObject(groups[1], FaceID)
	geompy.addToStudyInFather(domainCompound,groups[1], "Group interfaceface")

	# walls group
	#groups.append(geompy.CreateGroup(domainCompound, geompy.ShapeType["FACE"]))
	#FaceID = geompy.GetSubShapeID(domainCompound, faces[4])
	#geompy.AddObject(groups[2], FaceID)
	#geompy.addToStudyInFather(domainCompound,groups[2], "Group walls")

	# space group
	groups.append(geompy.CreateGroup(domainCompound, geompy.ShapeType["FACE"]))
	FaceID = geompy.GetSubShapeID(domainCompound, faces[2])
	geompy.AddObject(groups[2], FaceID)
	#FaceID = geompy.GetSubShapeID(domainCompound, faces[3])
	#geompy.AddObject(groups[2], FaceID)
	#FaceID = geompy.GetSubShapeID(domainCompound, faces[4])
	#geompy.AddObject(groups[2], FaceID)
	#FaceID = geompy.GetSubShapeID(domainCompound, faces[5])
	#geompy.AddObject(groups[2], FaceID)
	#FaceID = geompy.GetSubShapeID(domainCompound, faces[6])
	#geompy.AddObject(groups[2], FaceID)
	#FaceID = geompy.GetSubShapeID(domainCompound, faces[7])
	#geompy.AddObject(groups[2], FaceID)
	#geompy.addToStudyInFather(domainCompound,groups[2], "Group space")


	
	# ---------------------------------------------------------------- #
	# Mesh
	# ---------------------------------------------------------------- #

	# Face properties
	faceElems = []			# Number of elements per face
	triangleAreas = []		# Average area of a triangle
	linearSizes = []		# Average linear size of a triangle
	
	faceElems.append( numFacesBubble )
	faceElems.append( numFacesLiquidIO )
	#faceElems.append( numFacesLiquidW )
	
	triangleAreas.append( bubbleRadius*bubbleRadius*4*3.1415/faceElems[0] )
	triangleAreas.append( 4/3*3.1415*rSphere*rSphere/(faceElems[1]) )
	#triangleAreas.append( 2*3.1415*rBox*hBox/(faceElems[2]) )
	
	linearSizes.append( sqrt(4.0/sqrt(3)*triangleAreas[0]) )
	linearSizes.append( sqrt(4.0/sqrt(3)*triangleAreas[1]) )
	#linearSizes.append( sqrt(4.0/sqrt(3)*triangleAreas[2]) )	

	
	# Shell properties
	shell0Elems = numElemsBubble
	shell0Volume = 4.0/3.0*3.1415*pow(bubbleRadius,3.0)
	tet0Volume = shell0Volume/shell0Elems
	tet0linSize = pow( tet0Volume*12.0/sqrt(2), 1.0/3.0 )
	
	shell3Elems = numElemsLiquid
	shell3Volume = 4.0/3.0*3.1415*(pow(rSphere,3.0)-pow(bubbleRadius,3.0))
	tet3Volume = shell3Volume/shell3Elems
	tet3linSize = pow( tet3Volume*12.0/sqrt(2), 1.0/3.0 )
	
	print "tet0linSize = {0}".format(tet0linSize)
	print "tet3linSize = {0}".format(tet3linSize)
	
	# -----------------------------------------------------------------#
	
	# Compute source mesh for faces, i.e. a mesh on the innermost face
	# that is projected onto the other faces
	
	faceSourceMesh = smesh.Mesh(faces[0],"Face Source Mesh")
	faceSourceAlgo = faceSourceMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
	
	faceSourceAlgoParams = faceSourceAlgo.Parameters()
	faceSourceAlgoParams.SetMaxSize( linearSizes[0] )
	faceSourceAlgoParams.SetSecondOrder( 0 )
	faceSourceAlgoParams.SetOptimize( 1 )
	faceSourceAlgoParams.SetFineness( 2 )
	faceSourceAlgoParams.SetMinSize( linearSizes[0]/10 )
	faceSourceAlgoParams.SetQuadAllowed( allowQuadrangles )
	
	faceSourceMesh.Compute()

	#
	# Define face meshes for overall mesh
	#
	
	allMesh = smesh.Mesh(domainCompound, "Domain Compound")
	
	faceMeshes = []
	
	faceMeshes.append( allMesh.Projection1D2D(geom=faces[0]) )
	faceMeshes[0].SourceFace(faces[0],faceSourceMesh,None,None,None,None)
	
	faceMeshes.append( allMesh.Projection1D2D(geom=faces[1]) )
	faceMeshes[1].SourceFace(faces[0],faceSourceMesh,None,None,None,None)

	# Box
	faceMeshes.append( allMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=faces[2]) )
	#faceMeshes.append( allMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=faces[3]) )
	#faceMeshes.append( allMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=faces[4]) )
	#faceMeshes.append( allMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=faces[5]) )
	#faceMeshes.append( allMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=faces[6]) )
	#faceMeshes.append( allMesh.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=faces[7]) )

	# Box mesh parameters
	spaceFaceAlgoParams = faceMeshes[2].Parameters()
	spaceFaceAlgoParams.SetMaxSize( linearSizes[1] )
	spaceFaceAlgoParams.SetSecondOrder( 0 )
	spaceFaceAlgoParams.SetOptimize( 1 )
	spaceFaceAlgoParams.SetFineness( 2 )
	spaceFaceAlgoParams.SetMinSize( linearSizes[1]/10 )
	spaceFaceAlgoParams.SetQuadAllowed( allowQuadrangles )


	# Define shell meshes for overall mesh
	#
	
	shellMeshes = []
	
	# Shell 0 - bubble
	shellMeshes.append( allMesh.Tetrahedron(geom=shells[0]) )
	params = shellMeshes[0].Parameters()
	params.SetMaxSize( tet0linSize )
	params.SetSecondOrder( 255 )
	params.SetOptimize( 1 )
	params.SetFineness( 2 )
	params.SetMinSize( tet0linSize/20 )
	params.SetSecondOrder( 0 )
	# Viscous layers
	thickness = 0.2*bubbleRadius
	numberOfLayers = 3
	stretchFactor = 1.5
	ignoreFaces = []
	layersHyp = shellMeshes[0].ViscousLayers(thickness,numberOfLayers,stretchFactor,ignoreFaces)

	
	# Shell 3 - liquid
	shellMeshes.append( allMesh.Tetrahedron(geom=shells[1]) )
	params = shellMeshes[1].Parameters()
	params.SetMaxSize( tet3linSize )
	params.SetSecondOrder( 255 )
	params.SetOptimize( 1 )
	params.SetFineness( 2 )
	params.SetMinSize( tet3linSize/20 )
	params.SetSecondOrder( 0 )
	# Viscous layers
	thickness = 0.25*bubbleRadius
	numberOfLayers = 3
	stretchFactor = 1.5
	ignoreFaces = []
	#ignoreFaces = [faces[2],faces[3],faces[4],faces[5],faces[6],faces[7]]
	layersHyp = shellMeshes[1].ViscousLayers(thickness,numberOfLayers,stretchFactor,ignoreFaces)


	#
	# Define groups of mesh faces
	#
	
	allMesh.GroupOnGeom(groups[0],"interfaceShadow",SMESH.FACE)
	allMesh.GroupOnGeom(groups[1],"interface",SMESH.FACE)
	#allMesh.GroupOnGeom(groups[2],"walls",SMESH.FACE)
	allMesh.GroupOnGeom(groups[2],"space",SMESH.FACE)


	#
	# Compute and export mesh
	#
	allMesh.Compute()
	allMesh.ExportUNV(savePathUNV)
	allMesh.ExportSTL(savePathSTL)


	return 0

retVal = makeBubble()
if retVal!=0:
	print "Stopped. (Code {0})".format(retVal)
else:
	print "Done."

