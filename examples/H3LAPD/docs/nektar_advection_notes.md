Number of elements:
$N_x=5 ; N_y=5 ; N_z=10$

Number of faces:
$$
\begin{align}
N_{f,ext} &= 2N_xN_y + 2N_xN_z + 2N_yN_z &=~250 \\
N_{f,int} &= N_xN_y(N_z-1)+ N_xN_z(N_y-1)+ N_yN_z(N_x-1) &=~625 \\
N_{f,tot} &= N_{f,ext} + N_{f,int} &=~875
\end{align}
$$

Number of trace points: $N_{tr} =N_{f,tot}(N_{Modes}+1)^{N_{dim}-1}  = 875*(2+1)^2=7875$
