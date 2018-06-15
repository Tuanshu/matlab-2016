clear all

t = tcpip('140.112.172.1',23,'timeout',1);
fopen(t)
%%
p=1;
RawData=0;
try
    while 1
        p=p+1;
        RawData(p)=fread(t,1,'uint8');
        disp(p)
    end
end
%% Check number of Row
Return_Index=find(RawData==13);
if ~isempty(Return_Index)
    Number_of_Row=length(Return_Index);
    TheLastRow=RawData(Return_Index(end):end);
    CharTheLastRow=char(TheLastRow);
    CharData=char(RawData);
    CharTheLastRow

else
    CharData=char(RawData);
    CharTheLastRow=CharData;
    Number_of_Row=1;
    CharData 
end

%% Write ID and PW

ID='attacksoil';
fwrite(t,[ID char(13)]);

PW='as0228'
fwrite(t,[PW char(13)])


fwrite(t,[char(13)])
%%
p=1;
RawData=0;
while RawData(max(p,1))~=8
    p=p+1;
    RawData(p)=fread(t,1,'uint8');
    disp(p)
end
%% Check number of Row
Return_Index=find(RawData==13);
Number_of_Row=length(Return_Index);
TheLastRow=RawData(Return_Index(end):end);
CharData=char(RawData);
CharTheLastRow=char(TheLastRow);
CharTheLastRow
%fclose(t)


% 
% fid=fopen('D:\161102_TCPIP test\test.txt','w+')
% fwrite(fid,QQQ_Char,'char')
% fclose(fid)
% %%
% fwrite(t,'A')
% fclose(t)
% 
% response = t.read();
% 
% % 
% p=1;
% q=1;
% while p<4000
%     if QQQ(p)<128
%         QQQ_New(q)=QQQ(p);
%         p=p+1;
%         q=q+1;
%     else
%         QQQ_New(q)=QQQ(p)*256+QQQ(p+1);
%         p=p+2;
%         q=q+1;
%     end
%     disp(p);
% end
