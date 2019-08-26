load('ref_ribeira1bbb_reg1.mat');
refl = myShrinkImageByFactorD(reflectances,6);
main_fly(refl, 'scene7');
load('ref_braga1bb_reg1.mat');
refl = myShrinkImageByFactorD(reflectances,6);
main_fly(refl, 'scene6');