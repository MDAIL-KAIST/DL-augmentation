function new_atomtype = Func_atomnumber2atomtype(atomtype)
% J.Lee, KAIST, 2020

% non atom ->0
% ordering numbering

atomnumber=unique(atomtype);
new_atomtype=atomtype;
numbering=1;

    for i=1:1:length(atomnumber)
        if atomnumber(i)>0
            new_atomtype(new_atomtype==atomnumber(i))=numbering;
            numbering=numbering+1;
        else
            new_atomtype(new_atomtype==atomnumber(i))=0;
        end
    end
    
end


