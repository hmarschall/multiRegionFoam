Info << "Creating multi-region system\n" << endl;

autoPtr<multiRegionSystem> multiRegion(new multiRegionSystem(runTime));

regionTypeList& regions = multiRegion().regions();

bool LTS = false;

//forAll(regions, regionI)
//{
//    LTS = LTS || fv::localEulerDdt::enabled(regions[regionI]);
//}
