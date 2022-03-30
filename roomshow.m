function h_fig=roomshow(colobj)

%h_fig=figure(fignum);hold on;

for i=1:length(colobj.obj)
    [~,patchObj] = show(colobj.obj{i});    
    patchObj.FaceColor = colobj.color{i};
    patchObj.EdgeColor = 'none';
end

view(-90,90);axis([-1 11 -1 11 0 4]);
grid on;

end