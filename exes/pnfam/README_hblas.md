# The idea behind hamiltonian_blas.f90

In subroutine "meanfield", what we want to calculate can be written in the following way
$$ 
h_{ij} = \sum_{\sigma\sigma'} \int d^3r\ \psi^*_i\left(\vec{r},\sigma\right) \hat{h}\left(\vec{r}, \sigma\sigma'\right) \psi_j\left(\vec{r},\sigma'\right)
$$
where $\sigma$ and $\sigma'$ are spins. Only one coordinate $\vec{r}$ is involved because of local density approximation. 

To obtain the equation that can be conveniently used in pnFAM, we define a 5-dimensional vector wave function as 
$$
\vec{\psi}_i = \left(\psi_i, \frac{d\psi_i}{dr}, nl(i)\cdot y \cdot \psi_i, \frac{d\psi_i}{dz}, \nabla^2 \psi_i\right)
$$
where the third component is related to the derivative over angle $\phi$. 
Correspondingly the matrix element is calculated as 
$$
h_{ij} = \sum_{\sigma\sigma'} \int d^3r\ \vec{\psi}^*_i\left(\vec{r},\sigma\right) \odot \overleftrightarrow{h}\left(\vec{r}, \sigma\sigma'\right) \odot \vec{\psi}_j\left(\vec{r},\sigma'\right)
$$
where $\odot$ is the dot product in the newly-defined vector space and $\overleftrightarrow{h}$ is now a rank-2 tensor. 

In the code $\overleftrightarrow{h}$ is first constructed from local densities (see tables at the end); then 
$$
\vec{\psi}_j^h\left(\vec{r},\sigma\right) = \sum_{\sigma'} \overleftrightarrow{h}\left(\vec{r}, \sigma\sigma'\right) \odot \vec{\psi}_j\left(\vec{r},\sigma'\right)
$$
is calculated, and finally 
$$
h_{ij} = \sum_{\sigma} \int d^3r\ \vec{\psi}^*_i \left(\vec{r},\sigma\right) \odot \vec{\psi}^h_j\left(\vec{r},\sigma\right)
$$
is computed via matrix-matrix multiplication (dgemm) where $(\vec{r},\sigma)$ is the index to sum over. By sorting the basis states in each block according to their spins, we can easily separate calculations with different $\sigma$ (this also has something to do with parameter "LDC" in "dgemm", which is discussed at [Purpose of LDA argument in BLAS dgemm?](https://stackoverflow.com/questions/8206563/purpose-of-lda-argument-in-blas-dgemm)). 

Pairing part and density matrix 
$$
\rho\left(\vec{r}\sigma;\vec{r}'\sigma'\right) = \sum_{ij} \psi_i\left(\vec{r}\sigma\right) \rho_{ij} \psi^*_j\left(\vec{r}'\sigma'\right)
$$
can be computed in a similar way. However, only local densites are calculated instead of the full density matrix in coordinate space. 

Above idea is based on my experience with Sky3D, a HF+BCS+TDHF solver in 3D Cartesian coordinate space. Later I notice that HFODD paper also uses similar notations but things are handled in a different way in HFODD. 

## Appendix: tables about how $\overleftrightarrow{h}$ is constructed
For all the tables below, 
d2_a (full) = d2_a - xla^2 \* y^2 \* wf_a, d2_b (full) = d2_b -xlb^2 \* y^2 \* wf_b; 
a factor of 2 should be multiplied in the end. 

### ns(a) == ns(b)

<table>
   <tr>
      <td></td>
      <td>wf_b</td>
      <td>dr_b</td>
      <td>y*xlb*wf_b</td>
      <td>dz_b</td>
      <td>d2_b (full)</td>
   </tr>
   <tr>
      <td rowspan="5">wf_a</td>
      <td>2*crho*rho</td>
      <td>-i*cj*jr</td>
      <td>cj*jp</td>
      <td>-i*cj*jz</td>
      <td>2*cdrho*rho</td>
   </tr>
   <tr>
      <td>2*cs*sz*ns</td>
      <td>-crdj*(tjpz-tjzp)</td>
      <td>+i*csdj*jr*ns</td>
      <td>-crdj*(tjrp-tjpr)</td>
      <td>2*cds*sz*ns</td>
   </tr>
   <tr>
      <td>ct*tz*ns</td>
      <td>-csdj*jp*ns</td>
      <td>-i*crdj*(tjzr-tjrz)</td>
      <td>2*cgs*gs*ns</td>
      <td></td>
   </tr>
   <tr>
      <td>cf*fz*ns</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>ctau*tau</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="4">dr_a</td>
      <td>+i*cj*jr</td>
      <td>4*cdrho*rho</td>
      <td>crdj*rho*ns</td>
      <td>+i*csdj*sp</td>
      <td></td>
   </tr>
   <tr>
      <td>-crdj*(tjpz-tjzp)</td>
      <td>4*cds*sz*ns</td>
      <td>csdj*sz</td>
      <td>cf*sr*ns/2</td>
      <td></td>
   </tr>
   <tr>
      <td>-csdj*jp*ns</td>
      <td>ctau*rho</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>ct*sz*ns</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="4">y*xla*wf_a</td>
      <td>cj*jp</td>
      <td>crdj*rho*ns</td>
      <td>4*cdrho*rho</td>
      <td>-csdj*sr</td>
      <td></td>
   </tr>
   <tr>
      <td>-i*csdj*jr*ns</td>
      <td>csdj*sz</td>
      <td>4*cds*sz*ns</td>
      <td>-i*cf*sp*ns/2</td>
      <td></td>
   </tr>
   <tr>
      <td>+i*crdj*(tjzr-tjrz)</td>
      <td></td>
      <td>ctau*rho</td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td></td>
      <td>ct*sz*ns</td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="4">dz_a</td>
      <td>+i*cj*jz</td>
      <td>-i*csdj*sp</td>
      <td>-csdj*sr</td>
      <td>4*cdrho*rho</td>
      <td></td>
   </tr>
   <tr>
      <td>-crdj*(tjrp-tjpr)</td>
      <td>cf*sr*ns/2</td>
      <td>+i*cf*sp*ns/2</td>
      <td>4*cds*sz*ns</td>
      <td></td>
   </tr>
   <tr>
      <td>2*cgs*gs*ns</td>
      <td></td>
      <td></td>
      <td>cf*sz*ns ctau*rho</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td></td>
      <td></td>
      <td>ct*sz*ns</td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="2">d2_a (full)</td>
      <td>2*cdrho*rho</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>2*cds*sz*ns</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
</table>

Tensor part: 
<table>
   <tr>
      <td></td>
      <td>wf_b</td>
      <td>dr_b</td>
      <td>y*xlb*wf_b</td>
      <td>dz_b</td>
   </tr>
   <tr>
      <td rowspan="4">wf_a</td>
      <td></td>
      <td>+i*ctj1*(tjzr-tjrz)*ns</td>
      <td>ctj1*(tjpz-tjzp)*ns</td>
      <td>-i*ctj0*(tjrr+tjpp+tjzz)*ns</td>
   </tr>
   <tr>
      <td></td>
      <td>-i*ctj2*(tjrz+tjzr)*ns/2</td>
      <td>ctj2*(tjpz+tjzp)*ns/2</td>
      <td>+i*ctj2*(2*tjrr-tjpp-tjzz)*ns/9</td>
   </tr>
   <tr>
      <td></td>
      <td></td>
      <td></td>
      <td>+i*ctj2*(2*tjpp-tjzz-tjrr)*ns/9</td>
   </tr>
   <tr>
      <td></td>
      <td></td>
      <td></td>
      <td>-2i*ctj2*(2*tjzz-tjrr-tjpp)*ns/9</td>
   </tr>
   <tr>
      <td rowspan="2">dr_a</td>
      <td>-i*ctj1*(tjzr-tjrz)*ns</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj2*(tjrz+tjzr)*ns/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="2">y*xla*wf_a</td>
      <td>ctj1*(tjpz-tjzp)*ns</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>ctj2*(tjpz+tjzp)*ns/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="4">dz_a</td>
      <td>+i*ctj0*(tjrr+tjpp+tjzz)*ns</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjrr-tjpp-tjzz)*ns/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjpp-tjzz-tjrr)*ns/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+2i*ctj2*(2*tjzz-tjrr-tjpp)*ns/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
</table>

### ns(a) == +1, ns(b) == -1
<table>
   <tr>
      <td></td>
      <td>wf_b</td>
      <td>dr_b</td>
      <td>y*xlb*wf_b</td>
      <td>dz_b</td>
      <td>d2_b (full)</td>
   </tr>
   <tr>
      <td rowspan="3">wf_a</td>
      <td>2*cs*(sr-i*sp)</td>
      <td>-i*csdj*jz</td>
      <td>-i*csdj*jz</td>
      <td>csdj*(jp+i*jr)</td>
      <td>2*cds*(sr-i*sp)</td>
   </tr>
   <tr>
      <td>ct*(tr-i*tp)</td>
      <td>2*cgs*gs</td>
      <td>2*cgs*gs</td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>cf*(fr-i*fp)</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="3">dr_a</td>
      <td>-i*csdj*jz</td>
      <td>cf*sr</td>
      <td>cf*(sr+i*sp)/2</td>
      <td>cf*sz/2</td>
      <td></td>
   </tr>
   <tr>
      <td>2*cgs*gs</td>
      <td>ct*(sr-i*sp)</td>
      <td></td>
      <td>crdj*rho</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>4*cds*(sr-i*sp)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="3">y*xla*wf_a</td>
      <td>+i*csdj*jz</td>
      <td>-cf*(sr+i*sp)/2</td>
      <td>-i*cf*sp</td>
      <td>-cf*sz/2</td>
      <td></td>
   </tr>
   <tr>
      <td>-2*cgs*gs</td>
      <td></td>
      <td>ct*(sr-i*sp)</td>
      <td>-crdj*rho</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td></td>
      <td>4*cds*(sr-i*sp)</td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="2">dz_a</td>
      <td>csdj*(jp+i*jr)</td>
      <td>cf*sz/2</td>
      <td>cf*sz/2</td>
      <td>ct*(sr-i*sp)</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>-crdj*rho</td>
      <td>-crdj*rho</td>
      <td>4*cds*(sr-i*sp)</td>
      <td></td>
   </tr>
   <tr>
      <td>d2_a (full)</td>
      <td>2*cds*(sr-i*sp)</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
</table>

Tensor part: 
<table>
   <tr>
      <td></td>
      <td>wf_b</td>
      <td>dr_b</td>
      <td>y*xlb*wf_b</td>
      <td>dz_b</td>
   </tr>
   <tr>
      <td rowspan="6">wf_a</td>
      <td></td>
      <td>-ctj1*(tjrp-tjpr)</td>
      <td>-ctj1*(tjrp-tjpr)</td>
      <td>-i*ctj1*(tjzr-tjrz)</td>
   </tr>
   <tr>
      <td></td>
      <td>-ctj2*(tjrp+tjpr)/2</td>
      <td>+ctj2*(tjrp+tjpr)/2</td>
      <td>-i*ctj2*(tjrz+tjzr)/2</td>
   </tr>
   <tr>
      <td></td>
      <td>-i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td>-i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td>+ctj1*(tjpz-tjzp)</td>
   </tr>
   <tr>
      <td></td>
      <td>-2i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td>+i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td>-ctj2*(tjpz+tjzp)/2</td>
   </tr>
   <tr>
      <td></td>
      <td>+i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td>-2i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>+i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td>+i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="6">dr_a</td>
      <td>+ctj1*(tjrp-tjpr)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+ctj2*(tjrp+tjpr)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+2i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="6">y*xla*wf_a</td>
      <td>-ctj1*(tjrp-tjpr)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+ctj2*(tjrp+tjpr)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-2i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="4">dz_a</td>
      <td>+i*ctj1*(tjzr-tjrz)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj2*(tjrz+tjzr)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-ctj1*(tjpz-tjzp)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+ctj2*(tjpz+tjzp)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
</table>

### ns(a) == -1, ns(b) == +1
<table>
   <tr>
      <td></td>
      <td>wf_b</td>
      <td>dr_b</td>
      <td>y*xlb*wf_b</td>
      <td>dz_b</td>
      <td>d2_b (full)</td>
   </tr>
   <tr>
      <td rowspan="3">wf_a</td>
      <td>2*cs*(sr+i*sp)</td>
      <td>+i*csdj*jz</td>
      <td>-i*csdj*jz</td>
      <td>csdj*(jp-i*jr)</td>
      <td>2*cds*(sr+i*sp)</td>
   </tr>
   <tr>
      <td>ct*(tr+i*tp)</td>
      <td>2*cgs*gs</td>
      <td>-2*cgs*gs</td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>cf*(fr+i*fp)</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="3">dr_a</td>
      <td>+i*csdj*jz</td>
      <td>cf*sr</td>
      <td>-cf*(sr-i*sp)/2</td>
      <td>cf*sz/2</td>
      <td></td>
   </tr>
   <tr>
      <td>2*cgs*gs</td>
      <td>ct*(sr+i*sp)</td>
      <td></td>
      <td>-crdj*rho</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>4*cds*(sr+i*sp)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="3">y*xla*wf_a</td>
      <td>+i*csdj*jz</td>
      <td>cf*(sr-i*sp)/2</td>
      <td>+i*cf*sp</td>
      <td>cf*sz/2</td>
      <td></td>
   </tr>
   <tr>
      <td>2*cgs*gs</td>
      <td></td>
      <td>ct*(sr+i*sp)</td>
      <td>-crdj*rho</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td></td>
      <td>4*cds*(sr+i*sp)</td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="2">dz_a</td>
      <td>csdj*(jp-i*jr)</td>
      <td>cf*sz/2</td>
      <td>-cf*sz/2</td>
      <td>ct*(sr+i*sp)</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>crdj*rho</td>
      <td>-crdj*rho</td>
      <td>4*cds*(sr+i*sp)</td>
      <td></td>
   </tr>
   <tr>
      <td>d2_a (full)</td>
      <td>2*cds*(sr+i*sp)</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
</table>

Tensor part: 
<table>
   <tr>
      <td></td>
      <td>wf_b</td>
      <td>dr_b</td>
      <td>y*xlb*wf_b</td>
      <td>dz_b</td>
   </tr>
   <tr>
      <td rowspan="6">wf_a</td>
      <td></td>
      <td>+ctj1*(tjrp-tjpr)</td>
      <td>-ctj1*(tjrp-tjpr)</td>
      <td>-i*ctj1*(tjzr-tjrz)</td>
   </tr>
   <tr>
      <td></td>
      <td>+ctj2*(tjrp+tjpr)/2</td>
      <td>+ctj2*(tjrp+tjpr)/2</td>
      <td>-i*ctj2*(tjrz+tjzr)/2</td>
   </tr>
   <tr>
      <td></td>
      <td>-i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td>+i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td>-ctj1*(tjpz-tjzp)</td>
   </tr>
   <tr>
      <td></td>
      <td>-2i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td>-i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td>+ctj2*(tjpz+tjzp)/2</td>
   </tr>
   <tr>
      <td></td>
      <td>+i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td>+2i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td>+i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td>-i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="6">dr_a</td>
      <td>-ctj1*(tjrp-tjpr)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-ctj2*(tjrp+tjpr)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+2i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="6">y*xla*wf_a</td>
      <td>-ctj1*(tjrp-tjpr)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+ctj2*(tjrp+tjpr)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj0*(tjrr+tjpp+tjzz)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjrr-tjpp-tjzz)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+2i*ctj2*(2*tjpp-tjzz-tjrr)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-i*ctj2*(2*tjzz-tjrr-tjpp)/9</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td rowspan="4">dz_a</td>
      <td>+i*ctj1*(tjzr-tjrz)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+i*ctj2*(tjrz+tjzr)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>+ctj1*(tjpz-tjzp)</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
   <tr>
      <td>-ctj2*(tjpz+tjzp)/2</td>
      <td></td>
      <td></td>
      <td></td>
   </tr>
</table>