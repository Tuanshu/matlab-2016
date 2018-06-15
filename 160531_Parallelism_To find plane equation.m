clear all

First_index=[6 1];

Second_index=[2 1];


Matrx_Cell=cell(6,3);

%Measurement 20th (and 19th)
Matrx_Cell{1,1}= [-346.509 -199.951 -608.729;%Plane 000
                             -336.764 -230.461 -608.731;
                             -285.772 -223.378 -608.723;
                             -286.178 -191.189 -608.721];

Matrx_Cell{2,1}= [-301.332 -198.898 -596.61;%Plane 2
                            -301.223 -207.523 -596.632;
                            -307.818 -218.104 -596.678;
                            -318.68 -221.663 -596.712;
                            -330.59 -213.983 -596.708;
                            -333.2 -205.764 -596.689;
                            -332.992 -200.043 -596.675];

Matrx_Cell{3,1}= [-316.544 -209.669 -552.962;%Plane 11
                    -320.324 -208.179 -552.983;
                    -320.023 -204.003 -552.983;
                    -313.767 -203.923 -552.949;
                    -313.281 -206.558 -552.945;
                    -314.952 -209.54 -552.953];

Matrx_Cell{4,1}= [-129.645 -201.432 -609.443;%Plane 000000
                -147.163 -112.223 -609.41;
                -330.987 -180.518 -609.466;
                -332.249 -80.581 -609.421];

Matrx_Cell{5,1}= [-244.718 -137.018 -524.592;%Plane 111
                -240.849 -140.745 -524.591;
                -237.058 -136.873 -524.622;
                -240.698 -133.09 -524.625;
                -243.832 -134.471 -524.607;
                -239.085 -140.316 -524.599];

Matrx_Cell{6,1}= [-308.389 -213.306 -595.466;%Plane 5
                -325.094 -213.518 -595.466;
                -321.198 -187.933 -595.453;
                -307.49 -187.76 -595.454;
                -308.946 -72.466 -595.403;
                -301.072 -72.368 -595.401;
                -301.186 -63.38 -595.399;
                -310.957 -63.504 -595.4;
                -147.347 -88.441 -595.385;
                -147.456 -79.851 -595.382;
                -133.93 -83.025 -595.379;
                -129.347 -91.328 -595.384;
                -138.888 -203.232 -595.424;
                -156.351 -203.452 -595.433;
                -160.619 -225.398 -595.437;
                -143.067 -232.07 -595.426];

                        
                        
%Measurement 22th & 23th
Matrx_Cell{1,2}= [  103.16 -258.372 -608.487;%Plane 00000
                              106.605 -336.608 -608.486;
                              -176.23 -318.158 -608.567;
                              -178.248 -272.344 -608.572];

Matrx_Cell{2,2}= [-99.029 -211.926 -565.283;%Plane 8
                          -87.368 -220.671 -565.282;
                          -91.225 -225.407 -565.282;
                          -99.155 -225.756 -565.284;
                          40.506 -215.541 -565.24;
                          27.628 -216.109 -565.244;
                          28.102 -226.879 -565.245;
                          43.697 -226.191 -565.24;
                          41.255 -364.662 -565.147;
                          41.737 -375.606 -565.129;
                          29.121 -376.161 -565.131;
                          28.581 -363.892 -565.154;
                          -87.727 -364.571 -565.191;
                          -87.085 -379.157 -565.164;
                          -97.333 -379.576 -565.166;
                          -97.942 -365.715 -565.19];

Matrx_Cell{3,2}= [-27.145 -278.153 -538.005;%Plane 11111
                              -30.579 -278.305 -537.99;
                              -30.316 -284.286 -538.007;
                              -25.92 -284.092 -538.027;
                              -26.154 -278.77 -538.012;
                              -32.656 -281.398 -537.989];

Matrx_Cell{4,2}= [-13.046 -282.996 -581.722%Plane 2A
                              -16.911 -270.362 -581.728
                              -30.574 -267.082 -581.716
                              -41.175 -272.743 -581.693
                              -43.946 -288.858 -581.68];

Matrx_Cell{5,2}= [-28.863 -277.777 -538.067%Plane 1111B
                              -25.736 -279.343 -538.072
                              -25.532 -283.152 -538.067
                              -27.879 -285.289 -538.06
                              -31.446 -283.862 -538.055
                              -32.358 -281.055 -538.057
                              -30.312 -278.356 -538.064];

Matrx_Cell{6,2}= [-13.365 -289.677 -581.72;%Plane 2B
                          -13.807 -279.629 -581.722;
                          -17.898 -269.53 -581.728;
                          -26.17 -267.589 -581.721;
                          -38.794 -269.904 -581.699;
                          -43.376 -277.333 -581.685;
                          -42.854 -289.176 -581.682];
                        
                        
                        
                        
%Measurement (24 &) 25th
Matrx_Cell{1,3}= [  208.832 133.113 -537.799;%Plane 11111
                              211.946 130.636 -537.802;
                              211.089 126.935 -537.824;
                              207.283 126.489 -537.836;
                              204.677 129.826 -537.827;
                              205.9 132.266 -537.812;
                              207.88 133.284 -537.801];

Matrx_Cell{2,3}= [  176.404 315.72 -608.327;%Plane 00000
                              194.641 317.862 -608.325;
                              228.304 -24.087 -608.325;
                              180.25 -29.731 -608.326];

Matrx_Cell{3,3}= [  187.035 196.658 -608.331;%Plane 00000B
                              163.526 215.985 -608.332;
                              131.616 194.17 -608.333;
                              133.598 177.297 -608.333];

Matrx_Cell{4,3}= [  144.848 185.337 -596.269;%Plane 2
                              154.695 175.401 -596.238;
                              176.919 190.6 -596.271;
                              163.715 206.588 -596.308;
                              150.62 203.213 -596.322];

Matrx_Cell{5,3}= [  158.896 194.21 -552.651;%Plane 11111B
                              163.445 192.456 -552.654;
                              163.862 188.907 -552.638;
                              161.679 187.238 -552.635;
                              157.561 188.457 -552.634;
                              157.129 192.14 -552.644;
                              158.531 194.126 -552.652];

Matrix_Element_First=Matrx_Cell{First_index(1),First_index(2)};
Matrix_Element_Second=Matrx_Cell{Second_index(1),Second_index(2)};

X_array_First=Matrix_Element_First(:,1)';
Y_array_First=Matrix_Element_First(:,2)';
Z_array_First=Matrix_Element_First(:,3)';   
 
X_array_Second=Matrix_Element_Second(:,1)';
Y_array_Second=Matrix_Element_Second(:,2)';
Z_array_Second=Matrix_Element_Second(:,3)';


for p=1:2
    if p == 1
        XY=[X_array_First;Y_array_First;ones([1 length(X_array_First)])]';
        Z=Z_array_First';%-min(Z_max_Pos);

    elseif p == 2        
        XY=[X_array_Second;Y_array_Second;ones([1 length(X_array_Second)])]';
        Z=Z_array_Second';%-min(Z_max_Pos);
    end


    beta_matrix=zeros(size(XY,2),size(XY,1));
    error_array=zeros(size(XY,1),1);
    for q=1:size(XY,1)
        XY_used=XY;
        Z_used=Z;
        XY_used(q,:)=[];
        Z_used(q)=[];
        beta_matrix(:,q)=XY_used\Z_used;
        error_array(q)=sum((Z_used-XY_used*beta_matrix(:,q)).^2);
    end
    [value index]=min(error_array);

    beta=beta_matrix(:,index);
    if p == 1
        a_First=beta(1);
        b_First=beta(2);
        c_First=beta(3);
    else
        a_Second=beta(1);
        b_Second=beta(2);
        c_Second=beta(3);
    end
end
%% ���

Residual_First=a_First*X_array_First+b_First*Y_array_First+c_First-Z_array_First;
Residual_Second=a_Second*X_array_Second+b_Second*Y_array_Second+c_Second-Z_array_Second;
%%

a_Diff=a_First-a_Second;
b_Diff=b_First-b_Second;
c_Diff=c_First-c_Second;





disp(index);
disp(sprintf('a_First=%3f\n',a_First))
disp(sprintf('b_First=%3f\n',b_First))


disp(sprintf('a_Second=%3f\n',a_Second))
disp(sprintf('b_Second=%3f\n',b_Second))
disp(sprintf('c_Second=%3f\n',c_Second))


disp(sprintf('a_Diff=%3f\n',a_Diff))
disp(sprintf('b_Diff=%3f\n',b_Diff))

disp(sprintf('a_Diff x 20 x 1000 =%3f micron\n',a_Diff*20*1000))
disp(sprintf('b_Diff x 20 x 1000 =%3f micron\n',b_Diff*20*1000))

disp(sprintf('Total_Diff=%3f micron',(abs(a_Diff)*20+abs(b_Diff)*20)*1000))
