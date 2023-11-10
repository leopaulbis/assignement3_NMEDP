clear all; close all;

length    = '5'; %'5','10'%modify the length
typeElem  = 'Q'; %'P','Q'
orderElem = '1'; %'1','2'
fineLevel = '0'; %'0','1','2','3'

switch typeElem
    case 'P'
        geomElemType = 'triangle';
        binElemType = 1;
    case 'Q'
        geomElemType = 'quadrilateral';
        binElemType = 0;
    otherwise
        error('Not existent element type')
end

fileName = ['meshes/' length 'zanja' typeElem orderElem '_' fineLevel '_readed.mat'];
disp(fileName);
load(fileName)

fprintf('\nDislpay order %s %s mesh composed by %s nodes in the impermeable soil.\n',...
    orderElem,geomElemType,int2str(size(Tb_bottom,1)+1))
fprintf('Artificial boundary at %s meters\n',length)

lineWidth = 3;
figure(1)
hold on
axis equal
T_boun = Tb_bottom(:,[1 end]);
X_boun = X(T_boun(:),:);
plot(X_boun(:,1),X_boun(:,2),'r-','LineWidth',lineWidth)%impermeable soil
T_boun = Tb_artificial(:,[1 end]);
X_boun = X(T_boun(:),:);
plot(X_boun(:,1),X_boun(:,2),'c-','LineWidth',lineWidth)%artificial boundary
T_boun = Tb_symmetry(:,[1 end]);
X_boun = X(T_boun(:),:);
plot(X_boun(:,1),X_boun(:,2),'b-','LineWidth',lineWidth)%symetric boundary
T_boun = Tb_h1(:,[1 end]);
X_boun = X(T_boun(:),:);
plot(X_boun(:,1),X_boun(:,2),'y-','LineWidth',lineWidth)%pizo 1
T_boun = Tb_h2(:,[1 end]);
X_boun = X(T_boun(:),:);
plot(X_boun(:,1),X_boun(:,2),'g-','LineWidth',lineWidth)%pizo 2
T_boun = Tb_wall(:,[1 end]);
X_boun = X(T_boun(:),:);
X_xmin = min(X_boun(:,1)); X_xmax = max(X_boun(:,1));
X_ymin = min(X_boun(:,2)); X_ymax = max(X_boun(:,2)); X_ymid = X(Tb_h2(1),2);
plot([X_xmin X_xmin X_xmax X_xmax],[X_ymid X_ymin X_ymin X_ymax],'m-','LineWidth',2)%wall
plot(X(:,1),X(:,2),'k*')%nodes
if(typeElem=='Q')
    numEdges = 4;
else
    numEdges = 3;
end
for ielem=1:size(T,1)
    for iedge=1:numEdges
        n1 = iedge;
        n2 = iedge+1;
        if(n2>numEdges)
            n2 = 1;
        end
        xn1 = X(T(ielem,n1),:);
        xn2 = X(T(ielem,n2),:);
        plot([xn1(1) xn2(1)],[xn1(2) xn2(2)],'k-','LineWidth',1)%plot element 
    end
end
legend('Impermeable soil','Artificial boundary','Symmetric',...
    'Fixed piezometric 1','Fixed piezometric 2','Wall','Mesh nodes','Elements')


ielem=1;
figure; hold on
for inode=1:size(T,2)
    xnode = X(T(ielem,inode),1);
    ynode = X(T(ielem,inode),2);
    plot(xnode,ynode,'*')
    text(xnode+0.01,ynode+0.01,int2str(inode))
end




