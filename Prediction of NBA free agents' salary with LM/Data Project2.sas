ods html close;
ods preferences;
ods html;

options ps=80 ls=80;
title 'nba player salary analysis and prediction';

data a;
length Names $20 POS $5 OldT $5 NewT $5;
infile 'd:\2012-2013-4.csv' firstobs=2 dlm=',' obs=128;
input Names$ Salary TransS POS OldT NewT TeamSa Age ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF;
run;

/* first step: print data */
proc print data = a;
run;


/* second step: descriptive statistics */
proc means data = a n mean median std min max;
run;


/* third step: correlation matrix */
proc corr data = a;
run;


/* fourth step: box-cox transformation */
proc transreg ss2; model boxcox(Salary)=identity(TeamSa ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF);
run;


/* fifth step: test for normality */
proc univariate data = a;
  var TransS;
run;


/* sixth step: full model and residual diagnosis */
ods graphics on;
   
  proc reg data = a;
    model TransS = TeamSa ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF/vif;
run;


/* seventh step: model selection */
title 'model selection';
proc reg data=a plots=(criteria cp);
model TransS = TeamSa ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF/ selection=rsquare cp;
run;


proc reg data=a plots=(criteria sbc);
model TransS = TeamSa ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF/ selection=forward details=summary;
model TransS = TeamSa ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF/ selection=backward details=summary;
model TransS = TeamSa ReAge GP MPG PTS FG threeP FT ORPG DRPG RPG RPM APG APM BLKPG BLKPM STPG STPM TOPG TOPM PFPG PFPM ASTvsTO BLKvsPF STvsTO STvsPF/ selection=stepwise details=summary;
run;



/* eighth step: selected models: stb vif influence Lowess */
title 'model1: TeamSa ReAge MPG PTS RPM';
proc reg data = a
plots(only label)=(RStudentByLeverage CooksD DFFITS DFBETAS ObservedByPredicted residuals(smooth));
model TransS = TeamSa ReAge MPG PTS RPM/stb vif;
model TransS = TeamSa ReAge MPG PTS RPM/p r influence;
run;


title 'model2: TransS = TeamSa ReAge MPG PTS FG RPM APG STPM PFPM ASTvsTO STvsTO';
proc reg data = a
plots(only label)=(RStudentByLeverage CooksD DFFITS DFBETAS ObservedByPredicted residuals(smooth));
model TransS = TeamSa ReAge MPG PTS FG RPM APG STPM PFPM ASTvsTO STvsTO/stb vif;
model TransS = TeamSa ReAge MPG PTS FG RPM APG STPM PFPM ASTvsTO STvsTO/p r influence;
run;



/* ninth step: confirm normality of residuals */
title 'confirm normality of residuals';
proc reg;
model TransS = TeamSa ReAge MPG PTS RPM;
output out=regout1 r=resmath p=predmath;
run;

proc capability data=regout1 normal;
var resmath;
run;


/* tenth step: scatterplot matrix */
title 'scatterplot matrix';
ods html;
proc sgscatter data = a;
matrix Salary TransS TeamSa Age ReAge MPG PTS RPM/ diagonal=(histogram normal);
quit;
ods html close;
ods graphics off;
