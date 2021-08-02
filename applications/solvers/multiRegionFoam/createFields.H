// create fuel cell system with mesh...

Info << "Creating multi-region system\n" << endl;

multiRegionSystem multiRegion(mesh);

regionTypeList& regions = multiRegion.regions();

bool LTS = false;

//forAll(regions, regionI)
//{
//    LTS = LTS || fv::localEulerDdt::enabled(regions[regionI]);
//}
