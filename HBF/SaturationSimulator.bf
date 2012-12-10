skipCodeSelectionStep = 1;
ExecuteAFile 			(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR + "chooseGeneticCode.def");
ApplyGeneticCodeTable  (0);

ExecuteAFile 			(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR + "DescriptiveStatistics.bf");
ExecuteAFile 			(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR + "LocalMGREV.bf");
ExecuteAFile 			(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR + "ReadDelimitedFiles.bf");

/* configuration parameters 
----------------------------------------------------------
*/
ACCEPT_ROOTED_TREES = 1;


iterates 		  = 1; 
treeLengthSteps   = {{0.001,0.01,0.1,1.0}};
treeString		  = "((Taxon1:1.0,Taxon2:1.0):10.0,(Taxon3:1.0,Taxon4:1.0):0.0);";
/*
treeString		  = "(((((((((((Taxon_18_Sampled_94:1.214929130195376,Taxon_3_Sampled_93:0.21492913019537596):4.237650001824818,Taxon_0_Sampled_100:11.452579132020194):4.805248489054126,Taxon_14_Sampled_96:12.25782762107432):4.853148872485367,Taxon_5_Sampled_79:0.11097649355968642):2.256123267286661,((Taxon_8_Sampled_91:2.611300510714516,Taxon_23_Sampled_89:0.611300510714516):0.07955749522476019,Taxon_20_Sampled_89:0.6908580059392762):11.676241754907071):1.2356155116266123,Taxon_19_Sampled_77:1.6027152724729596):9.11240269987525,Taxon_17_Sampled_69:2.71511797234821):4.768184094344264,Taxon_16_Sampled_63:1.4833020666924739):9.020346041936655,Taxon_21_Sampled_56:3.503648108629129):0.5564041475962611,Taxon_10_Sampled_64:12.06005225622539):51.941495882532465,((((((((Taxon_4_Sampled_95:20.50312877280772,Taxon_13_Sampled_76:1.5031287728077203):0.5255901025118135,Taxon_6_Sampled_74:0.028718875319533765):1.0740738838662054,Taxon_22_Sampled_75:2.102792759185739):2.308735426183908,Taxon_7_Sampled_74:3.4115281853696473):6.822094240269017,Taxon_1_Sampled_68:4.233622425638664):4.239123551695414,(Taxon_9_Sampled_89:10.101652694384605,Taxon_15_Sampled_79:0.10165269438460456):19.371093282949474):15.829127978429959,Taxon_11_Sampled_62:18.301873955764037):2.3821463427432477,((Taxon_12_Sampled_51:0.8020885598169727,Taxon_24_Sampled_53:2.8020885598169727):6.415530959141542,Taxon_2_Sampled_54:10.217619518958514):2.4664007795487706):41.31752784025057);";
*/
stationaryFreqs = {
	{    0.25,    0.25,    0.25}
	{    0.25,    0.25,    0.25}
	{    0.25,    0.25,    0.25}
	{    0.25,    0.25,    0.25}
	};
	
	
rootSeq 		= "";
/* set to an empty string to simulate from random frequency based sequences */

AC 				= 0.5;
AT 				= 0.5;
CG 				= 0.5;
CT 				= 1.0;
GT 				= 0.5;	

_maxAllowedRate	 = 10;
pathToFEL		 = "FELoutput.csv";
saveReplicatesTo = "simulated_data/replicate";

/* end configuration parameters 
--------------------------------------------------------------
*/


dNdS 		= 			((ReadCSVTable		(pathToFEL, 1))["1"])["1"];
sites 		= 			Rows(dNdS);	

if (Abs (rootSeq) && Abs (rootSeq) != sites*3)
{
	fprintf (stdout, "[ERROR: Incompatible lengths of rootSeq and the matrix of dN/dS estimates]");
}

dNdSMatrix	=			{sites,2};
for (k = 0; k < sites; k = k+1)
{
	dNdSMatrix[k][0] = Min (_maxAllowedRate, dNdS[k][1]);
	dNdSMatrix[k][1] = Min (_maxAllowedRate, dNdS[k][2]);
}

meanDS	= (({1,sites}["1"])*dNdSMatrix[-1][0])[0]/sites;
meanDN	= (({1,sites}["1"])*dNdSMatrix[-1][1])[0]/meanDS/sites;
dNdS	= dNdSMatrix * (1/meanDS);

UseModel		 	(USE_NO_MODEL);

Tree				ScalingTree = treeString;

PopulateModelMatrix ("MG94Q", stationaryFreqs);
codonFreqs			= BuildCodonFrequencies (stationaryFreqs);
Model MG94 			= (MG94Q, codonFreqs, 0);


GetString		(branchLengthExpression, MG94, -1);
GetInformation  (mInfo, MG94);

solveFor		= "scaler";
locCount		= Columns (mInfo);
/*fprintf			(stdout, "Defined a model on ", locCount, " local variables:\n");*/
for (k = 0; k < locCount; k = k+1)
{
	branchLengthExpression = branchLengthExpression ^ {{mInfo[k],mInfo[k]+"*"+solveFor}};
}

Tree			T = treeString;
treeBranchNames   = BranchName (T,-1);
branchCount		  = Columns (treeBranchNames) - 1;

global alpha = 0.5;
alpha:>0.01;alpha:<100;
category c = (4, EQUAL, MEAN, 
				GammaDist(_x_,alpha,alpha), 
				CGammaDist(_x_,alpha,alpha), 
				0 , 
				1e25,
				CGammaDist(_x_,alpha+1,alpha)
			 );
			 
			  	 
global			  ACM = AC;
global			  ATM = AC;
global			  CGM = AC;
global			  CTM = CT;
global			  GTM = AC;



GTR_gamma = {{*,ACM*mu*c,mu*c,ATM*mu*c}{ACM*mu*c,*,CGM*mu*c,CTM*mu*c}{mu*c,CGM*mu*c,*,GTM*mu*c}{ATM*mu*c,CTM*mu*c,GTM*mu*c,*}};
nucFreqs  = (stationaryFreqs * {{1}{1}{1}}) * (1/3);

Model GTR = (GTR_gamma, nucFreqs, 1);
	  	 
codons = {{"A","C","G","T"}
          {"3",GeneticCodeExclusions,"",""}};
               
	  	 
	  	 
for (branchLength = 0 ; branchLength < Columns (treeLengthSteps); branchLength = branchLength + 1)
{
	blDef = treeLengthSteps[branchLength];

	for (siteID = 0; siteID < sites; siteID = siteID + 1)
	{
		synRate 		= dNdS[siteID][0];
		nonSynRate		= dNdS[siteID][1];
		
		if (synRate == 0 && nonSynRate == 0)
		{
			setExpr = 0;
		}
		else
		{
			ExecuteCommands ("FindRoot (setExpr, " + branchLengthExpression+ "-"+ 3*blDef + "," + solveFor + ",0,10000)");
		}
		ExecuteCommands		(solveFor + "=" + setExpr);
		for 	(bID = 0; bID < branchCount; bID = bID + 1)
		{
			ExecuteCommands ("T." + treeBranchNames[bID] + ".synRate="+solveFor+"*"+dNdS[siteID][0]);
			ExecuteCommands ("T." + treeBranchNames[bID] + ".nonSynRate="+solveFor+"*"+dNdS[siteID][1]);
		}
		
        spawnFrom = iterates;       
        if (Abs(rootSeq))
        {
        	spawnFrom = "";
        	template = rootSeq [siteID*3][siteID*3+2];
        	
        	for (siteCopy = 0; siteCopy < iterates; siteCopy = siteCopy + 1)
        	{
        		spawnFrom * template;
        	}
        	spawnFrom * 0;
        }

		if (siteID == 0)
		{
			DataSet			 bigSet = Simulate (T,codonFreqs,codons,spawnFrom);
		}
		else
		{
			DataSet			 thisSet = Simulate (T,codonFreqs,codons,spawnFrom);		
			DataSet			 bigSet  = Concatenate (bigSet, thisSet);
		}
	}
	
	recoveredLengths = {iterates, branchCount};
	pairwiseIdentity = {iterates, 1};

	saveRep = saveReplicatesTo + "_" + blDef + "_raw";
	DataSetFilter thisReplicate = CreateFilter (bigSet,1);
	fprintf (saveRep, CLEAR_FILE, thisReplicate);


}
