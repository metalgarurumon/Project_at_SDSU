libname iq "~/827/Data/";

data data1;
set iq.b680proj;
run;

data data2;
set data1;
if mst=5 then never=1;
else never=0;
run;

proc logistic data=data2;
strata str;
class fndx (ref="0") chk (ref="1") never (ref="0")/param=ref;
model fndx=chk agmn wt never/influence;
output out=diagn h=leve reschi=pear predicted=fitt;
run;


data daign2;
set diagn;
npear=(fndx-fitt)/sqrt(fitt);
chisqu=npear**2;
betas=chisqu*(leve/(1-leve));
run;

proc gplot data=daign2;
plot chisqu*fitt;
plot betas*fitt;
plot leve*fitt;
run;

proc export data=daign2 outfile="~/827/Data/diag_680.csv" dbms=csv replace;
run;