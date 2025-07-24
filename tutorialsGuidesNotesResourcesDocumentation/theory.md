This is 
$$
\renewcommand{\u}[1]{\mathbf{#1}}
\newcommand{\dd}{\mathrm{d}}
$$

$$
M\ddot{u}(t)+ K u(t)= f(t)
$$
$$
\textsf{rhs}=f- K u
$$
```cpp
K.Mult(u, rhs)
add(F, -1.0, rhs, rhs); // un-comment this if prescribing non-zero tractions.
```

$$
\ddot{u}=M^{-1}\textsf{rhs}
$$
```cpp
cg.Mult(rhs, a); 
```

$$
\begin{align}
\frac{u_{n+1}-2u_n+u_{n-1}}{\Delta t^2}&=\ddot{u}_n,\\
u_{n+1}-2u_n+u_{n-1}&=\ddot{u}_n \Delta t^2,\\
u_{n+1}&=2u_n-u_{n-1}+\ddot{u}_n \Delta t^2,\\
\end{align}
$$
		
```
add(u, u, temp);										  
add(temp, -1.0, u_old, temp);							  
add(temp, brainSimProps.Δt * brainSimProps.Δt, a, u_new); 
```



$\epsilon= u/L$
$U_{elastic}=\frac{1}{2} \sigma\epsilon$ $A L$
$U_{elastic}=\frac{1}{2} E\epsilon^2$ $A L$
$U_{elastic}=\frac{1}{2} E A \frac{u^2}{L}$

$U_{s}=2A\gamma$
$U_{elastic}-U_{surface}=\frac{A}{2}\left({E L} \frac{u}{L^2}^2-4\gamma\right)$ 

$\frac{u}{L}=\sqrt{\frac{4\gamma}{E L}}$


