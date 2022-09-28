const REAL boxlength=1;
const int num_of_sites=1;
const int site_ids[2]=
{
1,
2
};
const int num_of_forms_per_site[2]=
{
1,
2
};
const int n_0_0 = 100;
const REAL l_0_0 = 1.0;
const REAL x_0_0[n_0_0] = 
{
0.139006129982,
0.932451488644,
0.0492772911777,
0.664039144643,
0.560151932725,
0.703371064995,
0.675557726642,
0.917637544205,
0.342242848619,
0.278454139594,
0.303845827074,
0.234756601804,
0.849280910321,
0.780230292803,
0.360448043812,
0.691423536765,
0.583179759962,
0.28910820372,
0.415976186135,
0.394143710672,
0.177908597824,
0.941186583358,
0.210282678281,
0.745453059408,
0.953464231602,
0.114838275892,
0.960630825887,
0.924590346783,
0.777912121678,
0.429262520319,
0.228683728842,
0.492050619341,
0.286327221854,
0.738777292393,
0.863244575403,
0.441593682693,
0.116342545275,
0.228553453199,
0.85763785738,
0.72205447296,
0.354273340459,
0.14738685231,
0.227327414934,
0.358267547026,
0.101872751665,
0.189272465697,
0.0677634274816,
0.0051418212668,
0.536035913284,
0.384474284186,
0.252208401393,
0.217223637454,
0.103925156202,
0.241069763394,
0.0307662332713,
0.0487913535175,
0.0640408995982,
0.162145170772,
0.418284898485,
0.185928558317,
0.152047712077,
0.800771667068,
0.733244251212,
0.785276297463,
0.229410521573,
0.297284578519,
0.616199305284,
0.472243323264,
0.186948910844,
0.128072401339,
0.766257600732,
0.519658775768,
0.831752054852,
0.51305877655,
0.466780407294,
0.619874987089,
0.0396741423173,
0.119829298395,
0.138364584006,
0.352183704972,
0.469831703418,
0.206686868257,
0.770692912374,
0.078025996018,
0.855199982633,
0.535987349185,
0.832568283726,
0.863226007082,
0.115397026163,
0.907057040777,
0.540588892004,
0.414621241766,
0.745193987475,
0.980616548469,
0.0840409502934,
0.848480477371,
0.483624978127,
0.92439874212,
0.088274938369,
0.710296939088
};
const REAL y_0_0[n_0_0] = 
{
0.67625005574,
0.804587900317,
0.227204613252,
0.904082203925,
0.688001658498,
0.614004344935,
0.615251418692,
0.0851989507103,
0.456512360377,
0.112235546541,
0.32710760587,
0.117326014397,
0.38305912661,
0.177290209918,
0.978271673757,
0.426565023633,
0.0925768357921,
0.205258967321,
0.152188884222,
0.920992463175,
0.375986944103,
0.636333946184,
0.0848361714112,
0.506040996975,
0.0370984761311,
0.76157813685,
0.5627215731,
0.76496801416,
0.746688402076,
0.698702415989,
0.254679293906,
0.253407853703,
0.397691512339,
0.986258674667,
0.142093557933,
0.00686777098724,
0.54504472981,
0.713644626821,
0.643031177305,
0.928841202329,
0.311859413306,
0.704698967922,
0.998920521618,
0.621393664828,
0.86023061538,
0.0587962011606,
0.635508032796,
0.454954545099,
0.94709871425,
0.599820644638,
0.7890386815,
0.0695372831402,
0.541458548796,
0.897110160015,
0.922636398492,
0.21354494558,
0.737040644945,
0.661414548019,
0.0917412323012,
0.131166049412,
0.944506945692,
0.36809899139,
0.190585659144,
0.0216128070278,
0.171144861186,
0.307507712888,
0.815152667462,
0.82825447095,
0.811771040155,
0.631848488554,
0.457357344143,
0.167737549873,
0.170539912591,
0.737908903528,
0.612432927443,
0.645307678624,
0.440393212136,
0.872225270422,
0.88772949374,
0.0211992174334,
0.382679210924,
0.086755719537,
0.354010874852,
0.411036383848,
0.915572574344,
0.255955137518,
0.330351487146,
0.078637447083,
0.577168973049,
0.559390877389,
0.904073181307,
0.0376253781357,
0.771071485977,
0.417519774181,
0.821456704451,
0.557190224343,
0.103834423211,
0.00232293655494,
0.188384359091,
0.299000563192
};
const REAL z_0_0[n_0_0] = 
{
0.481663699313,
0.20986255629,
0.785614005108,
0.106914056023,
0.147805571986,
0.750352378898,
0.605478326169,
0.975693035626,
0.32015426875,
0.408138455631,
0.533634364978,
0.458237871339,
0.589911578163,
0.72474388079,
0.00177040876504,
0.183901537971,
0.371021256973,
0.83267107065,
0.500064103271,
0.125164146208,
0.396464121699,
0.3134865122,
0.896315810815,
0.985074698894,
0.785484348599,
0.100838073571,
0.600337178656,
0.38227011232,
0.314981472815,
0.0807960402882,
0.063634413913,
0.364891580199,
0.257597138199,
0.724026983626,
0.890700939945,
0.888025581563,
0.445549152605,
0.532221640131,
0.673240504385,
0.229605023322,
0.26034269006,
0.0663685798131,
0.446236541954,
0.130921675239,
0.442422651104,
0.842601759379,
0.326145163872,
0.759928355172,
0.290504693428,
0.962402642855,
0.0785101213469,
0.121034477778,
0.577216200147,
0.969927929638,
0.129854270517,
0.709718959092,
0.666354051368,
0.867765221698,
0.20050958408,
0.287662207894,
0.859321839697,
0.0983872712626,
0.842832593287,
0.753190179929,
0.828297221718,
0.216592457688,
0.647708250732,
0.504752360048,
0.349406427272,
0.499215099532,
0.0858113812105,
0.25178273089,
0.642973316603,
0.151073442735,
0.300587069811,
0.32979196525,
0.105087405304,
0.983871264815,
0.143096124117,
0.974337349033,
0.774004613055,
0.269744276204,
0.797666656402,
0.740682417492,
0.989863699373,
0.0207548961452,
0.421854140502,
0.27608273815,
0.464818488799,
0.239399867468,
0.57180897528,
0.0308817382777,
0.186346327287,
0.590409408357,
0.86572912727,
0.627896197062,
0.24554404236,
0.987261285254,
0.920482058014,
0.295403853061
};
const REAL q_0_0[n_0_0] = 
{
-0.702916465067,
1.24147170397,
-1.30175247592,
2.042524707,
2.16237878365,
-0.0944655047537,
0.26287423844,
-1.5946766395,
0.723141425556,
-1.38256435033,
1.80638733402,
0.0862464425587,
2.30039749766,
0.829781296602,
2.95712046843,
1.43324507313,
0.307318823562,
-0.381988232941,
-0.192996789658,
-0.290670832426,
-2.09239320446,
1.27620141511,
1.21357339903,
-1.82187794093,
-0.0220725904899,
0.540984132292,
1.73449748426,
-1.85004509371,
-0.512069953633,
-2.25270115693,
1.02454232532,
-2.01511217459,
-1.0356323081,
0.911603904576,
-3.4952604336,
-0.586386625659,
-2.06366389124,
1.82096079149,
-0.424814766647,
-0.564986583213,
0.308987812473,
-2.1727886092,
2.48962576084,
-0.352138868136,
1.6174396357,
0.640678888779,
0.801417193404,
0.113646071954,
1.86075582707,
0.717212111569,
-1.07135217988,
0.312049092069,
-1.91541523858,
-0.922462349077,
1.03723384104,
1.75713794689,
0.967309032847,
-1.18847890007,
0.863316141074,
-1.53907374385,
2.261296345,
-0.746769955561,
1.99561190783,
-0.340861519824,
-1.48271580016,
-0.854904301396,
-1.44157741377,
-1.21334463218,
0.287952016145,
1.74216099563,
1.71819436158,
-1.21426867806,
-3.43854120468,
-1.26353009974,
-0.726407056685,
-0.887271254751,
1.64215745082,
-1.33381518766,
-1.3242259927,
-1.39225685482,
1.96208148694,
1.52092335118,
-0.455042506839,
1.63590573993,
0.361357750932,
-0.798875648806,
0.674964892862,
-0.81392140582,
-1.80565434624,
1.37831934287,
0.951567616641,
-1.16734428831,
1.27422890494,
-1.35625123102,
-2.63707571517,
-1.01897686518,
1.74893300842,
1.36683654618,
-1.33985598509,
-1.25816310222
};
const REAL **x_site0 = new const REAL *[1];
const REAL **y_site0 = new const REAL *[1];
const REAL **z_site0 = new const REAL *[1];
const REAL **q_site0 = new const REAL *[1];
const int n_1_0 = 108;
const REAL l_1_0 = 0.5;
const REAL x_1_0[n_1_0] = 
{
0.732827206967,
0.821816954075,
0.812026884944,
0.587367540655,
0.547740350974,
0.750501051178,
0.851199362352,
0.714771086208,
0.669570752375,
0.848680388645,
0.805140266021,
0.599982691495,
0.740494842941,
0.711628797274,
0.852554101725,
0.732264519076,
0.847459059541,
0.610008422503,
0.59390515213,
0.723561303055,
0.74323924367,
0.690307782489,
0.720003289587,
0.922801961573,
0.799658314062,
0.553609823499,
0.69412409522,
0.566044387324,
0.67644050936,
0.906566235132,
0.663948502178,
0.840756891553,
0.883855395137,
0.789612240426,
0.698492013835,
0.674204382176,
0.785633532071,
0.649343447856,
0.724583721641,
0.852879158818,
0.707036161001,
0.642145398392,
0.700512045107,
0.757304217105,
0.598826900752,
0.746045321216,
0.842833733991,
0.610494805894,
0.669955404583,
0.685605541172,
0.711153946314,
0.671243797686,
0.721002850411,
0.596542709883,
0.485523139228,
0.747142620726,
0.930682932724,
0.603202147975,
0.751489670639,
0.600338117494,
0.695839893614,
0.775412047124,
0.746363393868,
0.687624051857,
0.907903449264,
0.899147981647,
0.6564157293,
0.695426818121,
0.726963299613,
0.511535805873,
0.664217439886,
0.607229661817,
0.859262199093,
0.53442375126,
0.874106133904,
0.861969118746,
0.826709506927,
0.827442269947,
0.722967286625,
0.72078633749,
0.752027201786,
0.840575184873,
0.971460056567,
0.812243647337,
0.650512000247,
0.777938696652,
0.83418291267,
0.681780159282,
0.69965566535,
0.794843416455,
0.721899155675,
0.797052775889,
0.710726654751,
0.741913720151,
0.832848605945,
0.82366160813,
0.765112952643,
0.703659388155,
0.749857905456,
0.733121167941,
0.959253366661,
0.889263758239,
0.688535787236,
0.702022290629,
0.629235910278,
0.661143752596,
0.743264789932,
0.696176936628
};
const REAL y_1_0[n_1_0] = 
{
0.164187921734,
0.0743870272758,
0.0173867637495,
0.0379337090306,
0.117950606293,
0.0331769218358,
0.162644829335,
0.0780808296404,
0.0697805260652,
0.0313971177222,
0.0852695478956,
0.0364698077522,
0.0951079913128,
0.192406761339,
0.0475172709436,
0.00865159314672,
0.0679559271572,
0.272498358075,
0.0328686756801,
0.124643347248,
0.0411591393971,
0.00792735259521,
0.0652876935524,
0.0383602259853,
0.100746874665,
0.0213556010912,
0.105615958559,
0.0165396722597,
0.0471232712753,
0.0266892471186,
0.0825501338246,
0.173528672633,
0.120025330425,
0.0404401259046,
0.0103869552796,
0.210013098121,
0.108074336406,
0.181769439289,
0.00145168406703,
0.00420774662928,
0.0804146685775,
0.226843624603,
0.00147869928956,
0.00660580108919,
0.0653705065508,
0.0251009092991,
0.28669555248,
0.0933187011476,
0.0590373402121,
0.107851832644,
0.000245617067815,
0.0480032360303,
0.062506810251,
0.0447718490761,
0.232298329467,
0.012277314748,
0.114412274615,
0.128542588276,
0.0958437096336,
0.100516117107,
0.109085797287,
0.0573489155739,
0.152864037242,
0.0248376680633,
0.00230090649434,
0.0519019682098,
0.0179815457627,
0.24283587168,
0.00847889232243,
0.000417347395963,
0.0219559142025,
0.108733171667,
0.0748355019117,
0.0606565844376,
0.0602259248916,
0.0389222417401,
0.13546850793,
0.0448700491934,
0.338289064021,
0.109614459485,
0.0814321840337,
0.0188077587876,
0.059722138013,
0.180859128861,
0.288028055072,
0.116992945827,
0.0323154671983,
0.000251441196473,
0.0685269249797,
0.157560729916,
0.137259064313,
0.0118873596723,
0.117288290236,
0.0271501344482,
0.00411213291147,
0.0753513001078,
0.0636042475797,
0.182029018854,
0.107180244727,
0.0745582922698,
0.04576477153,
0.17019835989,
0.00736267772865,
0.0811592282986,
0.0337290228067,
0.0444121020051,
0.0410723515349,
0.200049378377
};
const REAL z_1_0[n_1_0] = 
{
0.908880094273,
0.80534223616,
0.850171244369,
0.757304140899,
0.775615043766,
0.822180801118,
0.859884490858,
0.849411935487,
0.695725149844,
0.979236577837,
0.899604519071,
0.970974948431,
0.812011833813,
0.788043648349,
0.920239841603,
0.978449280093,
0.829328337548,
0.951927602044,
0.854507408237,
0.96930691244,
0.945291778338,
0.838061526331,
0.886622433349,
0.922065453661,
0.923367607311,
0.797174723805,
0.886695395771,
0.971554308986,
0.97722835825,
0.893562643968,
0.792513373382,
0.876222954326,
0.823324476653,
0.668571999514,
0.849724066678,
0.896538624547,
0.79245956435,
0.817582922941,
0.966805032334,
0.841938608102,
0.956148436141,
0.75060074296,
0.964224571628,
0.832788449639,
0.935198935081,
0.961867629184,
0.836733031525,
0.977816108243,
0.938817143565,
0.804548178102,
0.899732785477,
0.971046478348,
0.7542045887,
0.856845586215,
0.827880328956,
0.842250424431,
0.920492628786,
0.835373815409,
0.964778488594,
0.858744088698,
0.942025176723,
0.941108972417,
0.995179822723,
0.892261010598,
0.973943709473,
0.918365098843,
0.782240877868,
0.906978159457,
0.996086237171,
0.751321871425,
0.931932063874,
0.857118763455,
0.93343724404,
0.77760896298,
0.989676482145,
0.87607736411,
0.936491515934,
0.913721553699,
0.95736524478,
0.90381858254,
0.839068574306,
0.965030727356,
0.983672253775,
0.991612029996,
0.84562908994,
0.899789161436,
0.934318706605,
0.98578726119,
0.858230452523,
0.683998761213,
0.759093750674,
0.968905723701,
0.865058978866,
0.981720190696,
0.883473175582,
0.989837575639,
0.931815511136,
0.954826648818,
0.955468919106,
0.824910192087,
0.7736518894,
0.948958022479,
0.943074681017,
0.894598947969,
0.778593095934,
0.82278089305,
0.73813067956,
0.901033470068
};
const REAL q_1_0[n_1_0] = 
{
-0.296191964284,
-0.299175632768,
-0.304795679341,
-0.301220376365,
0.308747279311,
-0.302644454664,
0.301130527687,
-0.295017342898,
-0.297541211935,
-0.301882742704,
0.295480908352,
-0.301038501126,
0.297513606322,
-0.298779060488,
-0.308637826943,
-0.307857725094,
0.299009989816,
-0.298053352418,
0.297281884208,
-0.300834810082,
-0.298106158353,
0.297076335589,
-0.297849380926,
-0.301023348489,
0.306478980879,
-0.295095483367,
-0.290349829294,
0.304413971977,
0.30117639632,
-0.30738795898,
0.308907324514,
0.302477304301,
0.310532777624,
-0.304111531612,
-0.302010480721,
0.302534702313,
0.300871592553,
0.300884098567,
0.305516219436,
0.304180250977,
0.296728660175,
-0.292876708269,
0.301331111974,
0.298847371069,
0.297780933387,
0.294580081257,
0.305496233063,
0.292171274476,
-0.300503895845,
0.30132959484,
0.296160109226,
0.307435752461,
0.295779719637,
0.30245764761,
-0.296088411979,
0.296212095807,
0.297226119911,
-0.297023752325,
-0.303252246933,
-0.300831196068,
0.300699902562,
-0.303893642372,
0.297167469483,
-0.290082082396,
0.305539244008,
0.30464154591,
0.300270897139,
-0.305691239564,
-0.3017981973,
0.299639859401,
0.29953344086,
-0.300174975452,
0.304267072112,
0.303825125508,
0.293467659238,
-0.297816225328,
0.304719952204,
-0.295762103956,
-0.303937728218,
0.295448655892,
-0.295302378865,
-0.30377082764,
0.302526567582,
0.299579266284,
0.300724353794,
0.295917227961,
0.295257151226,
0.301490727958,
0.299176158667,
0.301344525924,
0.299636432296,
-0.298076302321,
0.294997941972,
0.307503133856,
-0.295344284569,
-0.301484079004,
0.298942392654,
0.301936496577,
0.300101074921,
-0.29626590125,
-0.300481552624,
0.301372681722,
-0.308636748816,
0.292002336619,
-0.302495191593,
0.297124166047,
0.306964932353,
-0.295862554449
};
const int n_1_1 = 109;
const REAL l_1_1 = 0.5;
const REAL x_1_1[n_1_1] = 
{
0.709093765343,
0.487653694511,
0.672134583239,
0.672239092366,
0.637156884262,
0.747301411038,
0.56621859601,
0.781253743222,
0.763116922678,
0.831772567219,
0.756545336984,
0.681972922006,
0.73296301567,
0.731959666431,
0.662123221974,
0.692713042947,
0.610063792516,
0.648067769939,
0.72460553784,
0.69624126772,
0.746692029919,
0.618212101414,
0.675343837356,
0.716700111367,
0.743036075783,
0.678781556854,
0.744821722031,
0.730072898544,
0.694846463662,
0.697637206399,
0.670474900185,
0.758709769219,
0.709402455213,
0.680225019782,
0.648167513383,
0.742968563597,
0.743828051094,
0.688814923522,
0.522063896896,
0.570101755962,
0.824886244309,
0.805807811325,
0.684045649141,
0.716993134785,
0.706863794597,
0.641275234011,
0.852341434454,
0.704075282656,
0.784169968833,
0.789037686101,
0.749156016165,
0.550346663413,
0.753730370438,
0.686099789491,
0.783356081722,
0.54997620482,
0.677434236181,
0.767527963345,
0.581646810123,
0.828037257083,
0.848802902762,
0.686884989483,
0.719248213655,
0.627765704313,
0.858303752065,
0.691401866441,
0.805908734143,
0.577403513136,
0.727853249842,
0.725990660816,
0.683369328167,
0.816750296483,
0.653950755633,
0.943606240213,
0.836834269108,
0.749793462366,
0.755997348394,
0.813365317442,
0.541076972672,
0.689421692336,
0.738650985515,
0.770156977677,
0.749756390836,
0.796997829139,
0.831915996841,
0.819706781826,
0.818810248595,
0.754282596525,
0.723928611246,
0.665226219173,
0.707972494928,
0.821302110358,
0.706084227542,
0.670320273752,
0.848026110157,
0.605888158996,
0.832529211534,
0.880249367872,
0.779709831488,
0.764713304349,
0.733134486849,
0.6278927998,
0.70657137236,
0.733252063032,
0.647030249933,
0.702921494203,
0.712196467437,
0.669690979682,
0.870159554022
};
const REAL y_1_1[n_1_1] = 
{
0.338723985816,
0.0783461999482,
0.0310969164181,
0.0898835882797,
0.08177755211,
0.195652690757,
0.198463208382,
0.00109675151493,
0.0959767682941,
0.0528935291127,
0.21794455663,
0.0733828628729,
0.192310060205,
0.0913527488662,
0.190520233978,
0.122660121144,
0.0321459988735,
0.0993247838646,
0.12211812918,
0.056419963616,
0.0602096738495,
0.100556504781,
0.046388230022,
0.102329358789,
0.0254571074709,
0.0515194593194,
0.0997236754293,
0.0905216372757,
0.0196511726716,
0.214394379506,
0.0048255042746,
0.147161870995,
0.105979482392,
0.0560270812885,
0.100633889525,
0.11834000077,
0.100456840551,
0.0209778934424,
0.123134834232,
0.0877315969359,
0.0724942266912,
0.134892966882,
0.0104092421419,
0.0791189071738,
0.222309728089,
0.0548815186876,
0.212438596016,
0.0541795110552,
0.0967602245961,
0.0197119328844,
0.208508355775,
0.0402018441272,
0.160011193214,
0.143477567721,
0.0770883825236,
0.0617566716895,
0.106542259353,
0.0142184232511,
0.210581763841,
0.106265158351,
0.0432482743188,
0.0115636468626,
0.0919646454325,
0.17579591787,
0.164132940549,
0.00290953691859,
0.124520489633,
0.125037359733,
0.0511219507663,
0.0634467880976,
0.0671227435467,
0.0138437623817,
0.0182919419182,
0.0356368460312,
0.0115634926929,
0.0380851934102,
0.157125894864,
0.0857497385844,
0.0564918526467,
0.031388971105,
0.0928324851382,
0.0522550902373,
0.0362599943992,
0.104717577276,
0.111424262853,
0.0786168867179,
0.034184907308,
0.155176852591,
0.189284759031,
0.0519729461103,
0.151238163696,
0.0973607686453,
0.0933895568056,
0.0990556405462,
0.125863587927,
0.00836828645622,
0.129959552735,
0.0958990220452,
0.171923412627,
0.324975082472,
0.0769395843666,
0.0384921664013,
0.010309646749,
0.00858044862283,
0.290621530377,
0.12279819745,
0.0514632978628,
0.0440327230639,
0.0790120736917
};
const REAL z_1_1[n_1_1] = 
{
0.994548805758,
0.954624871481,
0.916475345993,
0.886496071005,
0.914654472123,
0.935497902539,
0.991360977573,
0.771237423409,
0.799679500468,
0.985203165795,
0.965081366603,
0.894069622735,
0.949228689755,
0.978622387713,
0.93214973454,
0.777652819963,
0.926811195002,
0.915933070523,
0.717548726444,
0.805020122669,
0.807815903073,
0.957796582857,
0.953101937131,
0.991238521344,
0.797719993907,
0.985081752466,
0.836745680365,
0.954881215911,
0.8188024567,
0.849129829115,
0.865183153051,
0.969598939148,
0.980609307282,
0.905380219969,
0.931166259111,
0.880252930688,
0.920234729153,
0.819343709982,
0.945511932921,
0.913407471913,
0.99286278533,
0.873787145101,
0.964928096893,
0.868706380891,
0.904474509934,
0.954088668966,
0.860751044372,
0.72298919244,
0.986366573666,
0.875901190944,
0.924617129344,
0.95003685467,
0.907901733439,
0.743954714769,
0.911760411243,
0.939990487906,
0.958274974011,
0.879932413359,
0.856939335226,
0.875314019327,
0.836468109072,
0.988340843554,
0.950225874967,
0.75669423602,
0.937525562559,
0.916620509486,
0.93092976863,
0.733013148401,
0.929307195364,
0.94587786825,
0.988555962482,
0.944499952011,
0.892563591815,
0.931032794349,
0.958093109143,
0.905304232061,
0.860404626085,
0.906566658113,
0.85849296332,
0.987367473744,
0.777000077473,
0.964208839179,
0.80455739775,
0.789765125496,
0.708623363901,
0.944958391259,
0.914821849127,
0.89763121718,
0.673704265719,
0.855329742867,
0.906808224859,
0.887035927975,
0.993918074098,
0.928847339744,
0.967697269036,
0.88388007982,
0.91897870246,
0.784712507017,
0.934846641092,
0.912963027916,
0.989129432625,
0.992595318846,
0.94669945871,
0.764219351493,
0.922124129711,
0.968143670711,
0.923132732958,
0.693198280901,
0.918682328483
};
const REAL q_1_1[n_1_1] = 
{
0.298995530716,
0.291552683836,
-0.307450131491,
-0.307652539734,
0.288690395005,
-0.303702081141,
-0.304124911945,
-0.302850098944,
0.302943475675,
0.297358912897,
-0.304981281842,
0.299862465275,
0.302416247777,
-0.292103083586,
-0.298281455365,
0.288008751766,
0.299687983299,
0.291624906186,
0.29732607162,
0.299801298642,
0.299805758089,
-0.303597998329,
0.298099927994,
-0.305008755032,
0.303549997902,
-0.30514409557,
-0.305104931423,
-0.304065712028,
0.299423918287,
-0.307942145978,
-0.297467524787,
0.296685798156,
-0.307010251915,
0.293398737172,
-0.29838674484,
-0.298665748526,
0.296527329114,
0.291656303857,
0.306302792813,
0.29342824299,
0.297314838595,
0.306785223486,
-0.298288494912,
0.297480916955,
0.299005598819,
0.302961663575,
0.305057925842,
0.305309770905,
0.300774229722,
-0.295966640697,
-0.291786594011,
-0.299589517249,
0.303137839863,
0.304119803134,
0.303020659504,
-0.297825198604,
0.304196737579,
-0.305958914482,
-0.295130683752,
0.295707115714,
-0.299924980507,
-0.301803730278,
-0.292952344905,
-0.293559294689,
-0.291050062211,
-0.300253571898,
-0.297523981182,
0.294370127453,
-0.29965810236,
-0.29387987595,
0.297736083091,
0.304216216654,
0.299588239078,
-0.30923067694,
0.300927668636,
0.298118944926,
-0.30182969285,
0.296698858366,
0.302847210949,
-0.300627136419,
-0.308264624542,
-0.29412265959,
-0.306677054851,
0.295112858228,
0.303124606182,
0.294791721597,
0.29264205741,
0.294876370856,
-0.3006405896,
-0.301847144772,
-0.304884929075,
0.303683398755,
0.298863317675,
-0.304830479368,
-0.301301602424,
-0.30345300945,
0.296351675417,
-0.301267846797,
0.29279464773,
0.301089936551,
0.299031333393,
-0.300101561421,
-0.292414641257,
-0.300155720602,
-0.303927123019,
0.295275775408,
-0.300562874051,
-0.302485399871,
0.298135339869
};
const REAL **x_site1 = new const REAL *[2];
const REAL **y_site1 = new const REAL *[2];
const REAL **z_site1 = new const REAL *[2];
const REAL **q_site1 = new const REAL *[2];

const REAL ***x_data = new const REAL **[1];
const REAL ***y_data = new const REAL **[1];
const REAL ***z_data = new const REAL **[1];
const REAL ***q_data = new const REAL **[1];
int **n_data = new int *[1];
int *n_site0 = new int[1];
int *n_site1 = new int[2];

REAL **l_data = new REAL *[1];
REAL *l_site0 = new REAL[1];
REAL *l_site1 = new REAL[2];


int acu_data(){

x_site0[0] = x_0_0;
y_site0[0] = y_0_0;
z_site0[0] = z_0_0;
q_site0[0] = q_0_0;
n_site0[0] = n_0_0;
l_site0[0] = l_0_0;
x_site1[0] = x_1_0;
y_site1[0] = y_1_0;
z_site1[0] = z_1_0;
q_site1[0] = q_1_0;
n_site1[0] = n_1_0;
l_site1[0] = l_1_0;
x_site1[1] = x_1_1;
y_site1[1] = y_1_1;
z_site1[1] = z_1_1;
q_site1[1] = q_1_1;
n_site1[1] = n_1_1;
l_site1[1] = l_1_1;
x_data[0] = x_site0;
y_data[0] = y_site0;
z_data[0] = z_site0;
q_data[0] = q_site0;
n_data[0] = n_site0;
l_data[0] = l_site0;
x_data[1] = x_site1;
y_data[1] = y_site1;
z_data[1] = z_site1;
q_data[1] = q_site1;
n_data[1] = n_site1;
l_data[1] = l_site1;

return 1;
}
