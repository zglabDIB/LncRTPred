# LncRTPred
LncRTPred is a target prediction tool to predict the interaction between lncRNA and mRNA for both Human and Mouse species.
Webserver version corresponding to the LncRTPred tool can be accessed from http://bicresources.jcbose.ac.in/zhumur/lncrtpred/index.html

### Pre-requisite:
* The input dataset is provided in terms of sequences given in __fasta__ format corresponding to both lncRNA query and mRNA target for both species. The example of preformatted input sequence files (for both query and target) is provided in __Example__ folder.
* Sequence characters should be none other than __A, C, G, T__ or __U__.
* Input fasta format extension should be given as __.fa__; <br/>
  __.fasta__ format won’t be accepted.
* Input should be given in terms of __single-line fasta__ format corresponding to both lncRNA query and mRNA target files. Result from __multi-line fasta__ may contain error. 

__Software:__ For executing the codes, the required list of Softwares are given below.<br/>

__Python programming Environment:__
*	Anaconda software packages for python v3.8
*	pandas (min version: 1.1.3;  max version: 1.5.3; Version equal to or more than __2.0__ will fail to execute code.)
*	numpy  (min version: 1.19.1;  max version: 1.20.3.)
*	biopython version 1.78
*	scikit-Learn version 0.23.1  (min version: 0.23.1;  max version: 0.24.2; Version equal to or more than __1.0__ will fail to execute code.)
*	pickleshare version 0.7.5
*	lightgbm version 2.3.1

### Code Execution Procedure:
Note: After downloading from github, user may delete the __placeholder.txt__ file from both __Input_Data__ and __Output_Data__ folder in both Human and Mouse Directory before executing it.

__Human:__
* Place both the fasta files named __query_lncrna.fa__ and __target_mrna.fa__ in __Input_Data__ folder and execute it by the following code:
    * <i>__python predict_lncrna_mrna_human.py__</i>
*	After successful execution it will generate __Predicted_Val.csv__ in __Output_Data__ folder which provides the output in terms of predicted probabilities (in percentage) of successful interaction between lncRNA(s) and mRNA(s) given by __Decision Tree, k-nearest-neighbours, Random Forest__ and __LightGBM__ models.

___Mouse:___
* Place both the fasta files named __query_lncrna.fa__ and __target_mrna.fa__ in __Input_Data__ folder and execute it by the following code:
    * <i>__python predict_lncrna_mrna_mouse.py__</i>
* After the successful execution it will generate __Predicted_Val.csv__ in __Output_Data__ folder which provides the output in terms of predicted probabilities percentage of having successful interaction between lncRNA and mRNA given by __Decision Tree, k-nearest-neighbours, Random Forest__ and __LightGBM__ models.

### Final interaction score: 
* By default, LncRNA-mRNA interaction will be considered positive if prediction probability (in percentage) > 50 of at least 3 among all the 4 models, else negative.The last column named __Cumulative_Model_Score__ of the output file __Predicted_Val.csv__ file denotes the final interaction probability (in percentage). But, user may change the criterion according to their requirement.
