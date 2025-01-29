## Resubmission
This is a re-submission. In this version I have:


  
* provided an error message to provide guidance to the user if the BRW (2012)  procedure becomes highly collinear 
* replaced the use of the *index* function to eliminate the package dependence with the *zoo* package
* created the *FEVDandGFEVD* function that encompasses the functionalities of the previously available
  *FEVDsep*, *FEVDjoint*, *FEVDjointOrthoJLL*, *GFEVDsep*, *GFEVDjoint* and *GFEVDjointOrthoJLL*  functions.
* created the *OutputConstruction* function that encompasses the functionalities of the previously available *OutputConstructionSep*, *OutputConstructionJoint*.
* created the *TermPremiaDecomp* function that encompasses the functionalities of the previously available
  *TermPremiaDecompsep* and *TermPremiaDecompjoint*  functions.
* created the *FEVDandGFEVDgraphs* function that encompasses the functionalities of the
previously available *FEVDgraphsSep*, *FEVDgraphsJoint*, *FEVDgraphsJLLOrtho*,  *GFEVDgraphsSep*, *GFEVDgraphsJoint* and *GFEVDgraphsJLLOrtho* functions.
* created the *TPDecompGraph* function that encompasses the  functionalities of the previously available
  *TPDecompGraphsep* and *TPDecompGraphjoint*  functions.
* created the *FEVDandGFEVDgraphs*  function that encompasses the functionalities of the
  previously available *FEVDgraphsSep*,  *GFEVDgraphsSep*, *FEVDgraphsJoint*,  *GFEVDgraphsJoint*, *FEVDgraphsJLLOrtho* and *GFEVDgraphsJLLOrtho* functions.
* reshaped the *NumOutputs_Bootstrap* to accomodate the  functionalities of the previously available
  *OutputConstructionSep_BS* and  *OutputConstructionJoint_BS* functions.
* reshaped the *IRFandGIRFbs* to accomodate the functionalities of the previously available *IRFandGIRFbs_sepQ*,
  *IRFandGIRFbs_jointQ*, and *IRFandGIRFbs__jointQ_Ortho* functions.
* reshaped  the *FEVDandGFEVDbs* to accomodate the functionalities of the previously available
  *FEVDandGFEVDbs_sepQ*, *FEVDandGFEVDbs_jointQ*, and  *FEVDandGFEVDbs__jointQ_Ortho* functions.
