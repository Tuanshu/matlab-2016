clear all;

%題目為ax+by=z, a為2, b為1, 故有組合: X=1, 2 , 3,4   Y=3, 5, 7,9   Z=5, 9, 14
%(應為13),17

X=[1 3;2 5;3 7;4 9];
Y=[5;9;14;17];

beta_matrix=zeros(size(X,2),size(X,1));
error_array=zeros(size(X,1),1);
for p=1:size(X,1)
    X_used=X;
    Y_used=Y;
    X_used(p,:)=[];
    Y_used(p)=[];
    beta_matrix(:,p)=X_used\Y_used;
    error_array(p)=sum((Y_used-X_used*beta_matrix(:,p)).^2);
end
[value index]=min(error_array);

beta=beta_matrix(:,index);