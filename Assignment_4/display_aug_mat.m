function display_aug_mat(Aau)
[n,m]=size(Aau);
fprintf("$$ \n \\begin{bmatrix} \n")
for i=1:n
    for j=1:n
        fprintf("%8.2f & ",Aau(i,j));
    end
    fprintf(": & %8.2f \\\\ \n",Aau(i,n+1));
end



fprintf("\\end{bmatrix} \n $$ \n")
end

