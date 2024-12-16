## Resubmission
This is a re-submission. In this version, I have: 

* replaced the use of the *seqi* function to eliminate the package dependence with the *wrapr* package
* rewritten the overall code to allow to switch the *hablar* package from the *Imports* to the *Suggests* dependence section
* greatly simplified the previously named *getx* function and its *true2aux* sub-function (now named *Build_xvec* and  *GetAuxPara*, respectively)  
* greatly simplified the previously named *f_with_vectorized_parameters* function and its *update_para* and *aux2true* sub-functions (now named *Functionf_vectorized*, *Update_ParaList* , and  *GetTruePara*, respectively) 
* created the *VarianceExplained* function that encompasses the functionalities of the previously available *VarianceExplainedSep* and *VarianceExplainedJoint* functions
* created the *YieldsFit* function that encompasses the functionalities of the previously available *YieldsFitSep* and *YieldsFitJoint* functions
* created the *IRFandGIRF* function that encompasses the functionalities of the previously available *IRFsep*, *IRFjoint*, *IRFjointOrthoJLL*, *GIRFsep*, *GIRFjoint* and *GIRFjointOrthoJLL* functions.
* created the *IRFandGIRF_BS* function that encompasses the functionalities of the previously available *IRFsep_BS*, *IRFjoint_BS*, *IRFjointOrthoJLL_BS*, *GIRFsep_BS*, *GIRFjoint_BS* and *GIRFjointOrthoJLL_BS* functions.
* created the *IRFandGIRFgraphs* function that encompasses the functionalities of the previously available *IRFgraphsSep*, *IRFgraphsJoint*, *IRFgraphsJLLOrtho*, *GIRFgraphsSep*, *GIRFgraphsJoint* and *GIRFgraphsJLLOrtho* functions.
