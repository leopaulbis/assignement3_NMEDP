function [fileName]=mesh_pyzo(length,typeElem,orderElem,fineLevel)

%length    = '5'; %'5','10'%modify the length
%typeElem  = 'Q'; %'P','Q'
%orderElem = '2'; %'1','2'
%fineLevel = '2'; %'0','1','2','3'
orderElem=num2str(orderElem);
length=num2str(length);
fineLevel=num2str(fineLevel);
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

end
