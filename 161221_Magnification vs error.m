clear all

fobj=9;
fproj=50;
d=50;
D=0;
d_error=[50];

    for p=1:length(d_error)
    Image_Size=648*7.4; %micron

    Light_Image=[0;1];

    M_d=[1 d; 0 1];
    M_obj=[1 0; -1/fobj 1];
    M_proj=[1 0; -1/fproj 1];
    M_D=[1 D; 0 1];
    M_d_plus=[1 d+d_error(p)/2; 0 1];
    M_d_minus=[1 d-d_error(p)/2; 0 1];

    Light_Obj=M_obj*M_D*M_proj*M_d*Light_Image;
    Light_Obj_error_plus=M_obj*M_D*M_proj*M_d_plus*Light_Image;
    Light_Obj_error_minus=M_obj*M_D*M_proj*M_d_minus*Light_Image;

    Obj_Pos=Light_Obj(1)/Light_Obj(2)
    Obj_Pos_error_plus=Light_Obj_error_plus(1)/Light_Obj_error_plus(2)
    Obj_Pos_error_minus=Light_Obj_error_minus(1)/Light_Obj_error_minus(2)


    M=abs(Light_Obj(2)/Light_Image(2));

    M_error_plus=abs(Light_Obj_error_plus(2)/Light_Image(2));
    M_error_minus=abs(Light_Obj_error_minus(2)/Light_Image(2));

    Ratio=abs(M_error_plus-M_error_minus)/M;

    Size_Error(p)=abs(Image_Size/M_error_plus-Image_Size/M_error_minus);
    end

    plot(d_error,Size_Error);
    xlabel('CCD tolerance (mm)');
    ylabel('Object Space Size Error (micron)');