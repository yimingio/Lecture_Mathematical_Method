function display_matrix(A);
[n,m]=size(A);
fprintf("$$ \n \\begin{bmatrix} \n")
for i=1:n
    for j=1:m-1
        fprintf("%8.2f & ",A(i,j));
    end

    fprintf("%8.2f \\\\ \n",A(i,m));
end



fprintf("\\end{bmatrix} \n $$ \n")