%test linear solver

A = [ intval('[4,5]') intval('[-1,1]') intval('[1.5,2.5]') ;
         intval('[-0.5,0.5]') intval('[-7,-5]') intval('[1,2]') ;
          intval('[-1.5,-0.5]') intval('[-0.7,-0.5]') intval('[2,3]') ];
b =  [ intval('[3,4]'); intval('[0,2]'); intval('[3,4]') ];

X = verifylss(A,b)

pause

%test non-linear solver
xs = [0.05; 0.03]
[X , xs, k]= verifynlss(@ftest, xs , 'g', 'see')

pause 

%test gaussSeidelIteration
p =  [ intval('[-10,10]'); intval('[-10,10]'); intval('[-10,10]') ];

for i = 1:100
    X =  gaussSeidelIteration( A, b,  p);
    if ((i==1)||(i==2)||(i==5)||(i==10)||(i==20)||(i==100))
        X
    end
    p = X;
end


