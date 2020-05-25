function cur_Q = calcModularity(cur_w, cur_comm),

%Using definition for directed graphs.  From  V Nicosia et al J. Stat. Mech. (2009):
%  http://iopscience.iop.org/article/10.1088/1742-5468/2009/03/P03024/pdf

N = size(cur_w, 1);
m = sum(cur_w(:));
k_in = sum(cur_w, 1);
k_out = sum(cur_w, 2);

cur_Q = 0;
for i = 1:N,
    j_list = find(((cur_w(i, :) > 0)' | (cur_w(:, i) > 0)) & (cur_comm == cur_comm(i)));
    cur_Q = cur_Q + sum(cur_w(i, j_list) - k_in(i)*k_out(j_list)/m);
end
cur_Q = cur_Q/m;