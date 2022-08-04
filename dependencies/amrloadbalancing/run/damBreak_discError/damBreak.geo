//Variables

cellSize = 0.015;
convertToMeters = 0.146;


//pointBezeichener = newp; Point(pointBezeichner) = {x,y,z,cell size}
pointBL = newp; Point(pointBL) = {0*convertToMeters, 0*convertToMeters, 0*convertToMeters, cellSize};
pointBR = newp; Point(pointBR) = {4*convertToMeters, 0*convertToMeters, 0*convertToMeters, cellSize};
pointUL = newp; Point(pointUL) = {0*convertToMeters, 4*convertToMeters, 0*convertToMeters, cellSize};
pointUR = newp; Point(pointUR) = {4*convertToMeters, 4*convertToMeters, 0*convertToMeters, cellSize};

pointObstBL = newp; Point(pointObstBL) = {2*convertToMeters, 0*convertToMeters, 0*convertToMeters, cellSize};
pointObstBR = newp; Point(pointObstBR) = {2.16438*convertToMeters, 0*convertToMeters, 0*convertToMeters, cellSize};
pointObstUL = newp; Point(pointObstUL) = {2*convertToMeters, 0.32876*convertToMeters, 0*convertToMeters, cellSize};
pointObstUR = newp; Point(pointObstUR) = {2.16438*convertToMeters, 0.32876*convertToMeters, 0*convertToMeters, cellSize};

//lineBezeichner = newreg; Line(lineBezeichner) = {pointBezeichner1,pointBezeichner2};
lineL = newreg; Line(lineL) = {pointBL, pointUL};
lineR = newreg; Line(lineR) = {pointUR, pointBR};
lineU = newreg; Line(lineU) = {pointUL, pointUR};
lineBL = newreg; Line(lineBL) = {pointObstBL, pointBL};
lineBR = newreg; Line(lineBR) = {pointBR, pointObstBR};

lineObstL = newreg; Line(lineObstL) = {pointObstUL, pointObstBL};
lineObstR = newreg; Line(lineObstR) = {pointObstBR, pointObstUR};
lineObstU = newreg; Line(lineObstU) = {pointObstUR, pointObstUL};

//lineLoopBezeichner = newreg; Line Loop(lineLoopBezeichner) = {Liste mit lines};
lineLoop = newreg; Line Loop(lineLoop) = {lineL, lineU, lineR, lineBR, lineObstR, lineObstU, lineObstL, lineBL};

//planeSurfaceBezeichner = newreg; Plane Surface(planeSurfaceBezeichner) = lineLoopBezeichner;
planeSurface = newreg; Plane Surface(planeSurface) = lineLoop;

//extrude surface
Exlist[] = Extrude {0, 0, 0.01} {
		Surface{planeSurface};
		Layers{1};
		Recombine;
	} ;

Physical Volume("internalMesh") = {1};
Physical Surface("leftWall") = {23};
Physical Surface("atmosphere") = {27};
Physical Surface("rightWall") = {31};
Physical Surface("lowerWall") = {35, 39, 43, 47, 51};
Physical Surface("front") = {planeSurface};
Physical Surface("back") = {52};
