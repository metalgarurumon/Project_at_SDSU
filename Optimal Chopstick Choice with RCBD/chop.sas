
data a;

input efficiency trt block;

cards;
   19.55       1       1
   27.24       1       2
   28.76       1       3
   31.19       1       4
   21.91       1       5
   27.62       1       6
   29.46       1       7
   26.35       1       8
   26.69       1       9
   30.22       1      10
   27.81       1      11
   23.46       1      12
   23.64       1      13
   27.85       1      14
   20.62       1      15
   25.35       1      16
   28.00       1      17
   23.49       1      18
   27.77       1      19
   18.48       1      20
   23.01       1      21
   22.66       1      22
   23.24       1      23
   22.82       1      24
   17.94       1      25
   26.67       1      26
   28.98       1      27
   21.48       1      28
   14.47       1      29
   28.29       1      30
   27.97       1      31
   23.53       2       1
   26.39       2       2
   30.90       2       3
   26.05       2       4
   23.27       2       5
   29.17       2       6
   30.93       2       7
   17.55       2       8
   32.55       2       9
   28.87       2      10
   26.53       2      11
   25.26       2      12
   25.65       2      13
   29.39       2      14
   23.26       2      15
   24.77       2      16
   25.42       2      17
   23.65       2      18
   32.22       2      19
   18.86       2      20
   21.75       2      21
   23.07       2      22
   22.30       2      23
   27.04       2      24
   22.24       2      25
   24.87       2      26
   30.85       2      27
   21.15       2      28
   16.47       2      29
   29.05       2      30
   26.99       2      31
   21.34       3       1
   29.94       3       2
   32.95       3       3
   29.40       3       4
   22.32       3       5
   28.36       3       6
   28.49       3       7
   22.24       3       8
   36.15       3       9
   30.62       3      10
   26.53       3      11
   27.95       3      12
   31.49       3      13
   30.24       3      14
   24.80       3      15
   26.43       3      16
   29.35       3      17
   21.15       3      18
   29.18       3      19
   21.60       3      20
   25.39       3      21
   22.26       3      22
   24.85       3      23
   24.56       3      24
   16.35       3      25
   22.96       3      26
   25.82       3      27
   19.46       3      28
   23.60       3      29
   33.10       3      30
   27.13       3      31
   24.40       4       1
   25.88       4       2
   27.97       4       3
   24.54       4       4
   22.66       4       5
   28.94       4       6
   30.72       4       7
   16.70       4       8
   30.27       4       9
   26.29       4      10
   22.33       4      11
   24.85       4      12
   24.33       4      13
   24.50       4      14
   22.67       4      15
   22.28       4      16
   23.80       4      17
   25.36       4      18
   29.50       4      19
   20.19       4      20
   20.14       4      21
   21.09       4      22
   24.78       4      23
   24.74       4      24
   22.73       4      25
   21.08       4      26
   25.70       4      27
   19.79       4      28
   16.82       4      29
   31.15       4      30
   27.84       4      31
   22.50       5       1
   23.10       5       2
   28.26       5       3
   25.55       5       4
   16.71       5       5
   27.88       5       6
   31.07       5       7
   23.44       5       8
   28.82       5       9
   27.77       5      10
   24.54       5      11
   24.55       5      12
   27.78       5      13
   26.14       5      14
   23.44       5      15
   26.44       5      16
   27.47       5      17
   24.94       5      18
   29.68       5      19
   24.33       5      20
   25.42       5      21
   24.64       5      22
   22.78       5      23
   26.50       5      24
   18.71       5      25
   22.86       5      26
   25.09       5      27
   19.72       5      28
   17.05       5      29
   30.91       5      30
   25.92       5      31
   21.32       6       1
   26.18       6       2
   25.93       6       3
   28.61       6       4
   20.54       6       5
   26.44       6       6
   29.36       6       7
   19.77       6       8
   31.69       6       9
   24.64       6      10
   22.09       6      11
   23.42       6      12
   28.63       6      13
   26.30       6      14
   22.89       6      15
   22.68       6      16
   30.92       6      17
   20.74       6      18
   27.24       6      19
   17.12       6      20
   23.63       6      21
   20.91       6      22
   23.49       6      23
   24.86       6      24
   16.28       6      25
   21.52       6      26
   27.22       6      27
   17.41       6      28
   16.42       6      29
   28.22       6      30
   27.52       6      31

;
run;

run;




proc print data = a;

run;




proc means data = a maxdec = 3;

class trt;

var efficiency;

run;



proc glm data = a;

class block trt;

model efficiency = block;
random trt;

estimate 'short vs long ' trt 1 1 1 -1 -1 -1 / divisor = 3;

contrast 'short vs long ' trt 1 1 1 -1 -1 -1;

estimate '240mm vs 180mm ' trt -1 0 1 0 0 0;

estimate '240mm vs 210mm ' trt 0 -1 1 0 0 0;

estimate '240mm vs 270mm ' trt 0 0 1 -1 0 0;

estimate '240mm vs 300mm ' trt 0 0 1 0 -1 0;

estimate '240mm vs 330mm ' trt 0 0 1 0 0 -1;

contrast 'short vs long ' trt 1 1 1 -1 -1 -1;

run;

proc mixed data = a;

class block trt;
model efficiency = trt/solution;
random block; 
lsmeans trt/pdiff;

estimate 'short vs long ' trt 1 1 1 -1 -1 -1 / divisor = 3;
estimate '240mm vs 180mm ' trt -1 0 1 0 0 0;
estimate '240mm vs 210mm ' trt 0 -1 1 0 0 0;
estimate '240mm vs 270mm ' trt 0 0 1 -1 0 0;
estimate '240mm vs 300mm ' trt 0 0 1 0 -1 0;
estimate '240mm vs 330mm ' trt 0 0 1 0 0 -1;
run;
quit;

proc glmpower data=a;

class block trt;
model efficiency = trt block;
contrast 'short vs long ' trt 1 1 1 -1 -1 -1;

contrast '240mm vs 180mm ' trt -1 0 1 0 0 0;

contrast '240mm vs 210mm ' trt 0 -1 1 0 0 0;

contrast '240mm vs 270mm ' trt 0 0 1 -1 0 0;

contrast '240mm vs 300mm ' trt 0 0 1 0 -1 0;

contrast '240mm vs 330mm ' trt 0 0 1 0 0 -1;

power
stddev = 5
alpha = 0.025
ntotal = .
power = 0.8;
plot x=power min = .2 max =.9;

run;


Proc Univariate data=a Normal Plot ;
Var efficiency;
run;


proc boxplot data=a;
  plot efficiency*trt; 
      
run;


quit;


proc print data = diagnost;
run;

proc univariate data= diagnost freq normal plot normal noprint;
var resid;
title 'normal plot';
histogram/normal;
qqplot / normal;
run;

proc sgplot data = diagnost;
by trt;
scatter y=resid x=ybar;
title 'residual analyses';
run;
