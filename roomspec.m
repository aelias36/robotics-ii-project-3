function colobj = roomspec()

%
% definition file for objects in a scene
%

% constants
wlong=10;wshort=.5;wheight=4;
ex=[1;0;0];ey=[0;1;0];ez=[0;0;1];

% define for each object
k=1;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[wlong;wshort;wheight]; % size of object 
colobj.pos{k}=[wlong;-wshort;wheight]/2; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 0.5 0];
colobj.normaldir{k}=2;
% vector to surface:
% p0=colobj.pos{k}+colobj.ori{k}*(unitvec(colobj.normdir{k})*colobj.size{k}(colobj.normdir{k}))
% four corners:
% a=[1 2 3];a[colobj.normdir{k}]=[];
% b=[0;0;0];
% b(a(1))=+/-colobj.size{k}(a(1))/2;% b(a(2))=+/-colobj.size{k}(a(2))/2;
% corners=p0+colobj.ori{k}*b
% corner vectors: ecorner (+/-,+/1): (-/+,-/+) directions

% define for each object
k=2;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[wshort;wlong;wheight]; % size of object 
colobj.pos{k}=[-wshort;wlong;wheight]/2; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 0.5 0];

% define for each object
k=3;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[wlong;wshort;wheight]; % size of object 
colobj.pos{k}=[wlong/2;wlong+wshort/2;wheight/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 0.5 0];

% define for each object
k=4;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[wshort;wlong;wheight]; % size of object 
colobj.pos{k}=[wlong+wshort/2;wlong/2;wheight/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 0.5 0];

% define for each object
dL=1;dW=2;dH=.9;
desk1loc=[5;5];desk2loc=[2;8];

k=5;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[dL;dW;dH]; % size of object 
colobj.pos{k}=[desk1loc;;dH/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[.8 0 1];
% define for each object
k=6;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[dW;dL;dH]; % size of object 
colobj.pos{k}=[desk2loc;;dH/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[.8 0 1];

% define for each object
sR=.25;sL=.4;
k=7;
colobj.type(k)=2; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[sR;sL]; % size of object 
colobj.pos{k}=[4.2;4.2;sL/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[0 1 1];

% define for each object
k=8;
colobj.type(k)=2; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[sR;sL]; % size of object 
colobj.pos{k}=[4.2;5.8;sL/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[0 1 1];

% define for each object
k=9;
colobj.type(k)=2; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[sR;sL]; % size of object 
colobj.pos{k}=[2;7;sL/2]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[0 1 1];

% define for each object
k=10;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[1;2;2]; % size of object 
colobj.pos{k}=[9.5;8;1]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 .8 .8];

% define for each object
k=11;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[1;2;2]; % size of object 
colobj.pos{k}=[9.5;1;1]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 .8 .8];

% define for each object
k=12;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[2;1;2]; % size of object 
colobj.pos{k}=[4;0.5;1]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 .8 .8];

% define for each object
k=13;
colobj.type(k)=1; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[2;1.5;2]; % size of object 
colobj.pos{k}=[6;0.5;1]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[1 .8 .8];

% define for each object
k=14;
colobj.type(k)=2; % 1=Box, 2=Cylinder, 3=Sphere
colobj.size{k}=[.5,1.5]; % size of object 
colobj.pos{k}=[8;5;.75]; % location of box center
colobj.ori{k}=eye(3,3); % object orientation
colobj.color{k}=[.5 1 .5];

% 
for k=1:length(colobj.type)
    if colobj.type(k)==1 
            colobj.obj{k}=collisionBox(colobj.size{k}(1),...
            colobj.size{k}(2),colobj.size{k}(3));
    elseif colobj.type(k)==2
            colobj.obj{k}=collisionCylinder(colobj.size{k}(1),...
            colobj.size{k}(2));
    elseif colobj.type(k)==3
            colobj.obj{k}=collisionCylinder(colobj.size{k}(1));
    end
    colobj.obj{k}.Pose(1:3,4)=colobj.pos{k};
    colobj.obj{k}.Pose(1:3,1:3)=colobj.ori{k};
end

end