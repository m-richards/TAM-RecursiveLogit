function to_sparse_file(mat, fname)
[row, col, dat] = find(mat);
dlmwrite(fname,[row, col, dat], 'delimiter', '\t', 'newline', 'pc') 

end