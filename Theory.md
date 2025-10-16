# Constructing a helix using the _Joint Face Normal Method_

This notebook will implement the following algorithm:

Given :
- the length $L$ of an object,
- the angles between the axis and the normal vectors of the joint faces $\theta/2$, where the angles are assumed to be the same but in opposite directions and the two face normals and the axis are co-planar. Note that the total angle between the planes defined by the joint faces are $\theta$, which is consistent with the reference,
- the twist angle $\tau$. This is the same angle described in the reference.

All relevant properties of the matrix will be computed. Then later, practical methods to find the right parameters for the welding application.


## Homogeneous coordinates and rigid transformations

In **homogeneous coordinates** in 3D, a space is described by a vector $(x, y, z, w)$. The corresponding cartesian coordinates are $(x/w, y/w, z/w)$. The point $(x, y, z, 0)$ signifies the point 'at infinity, so meaning a direction.

**Rigid body transormations** are transormations, where all distances are preserved. They include rotations and translations, which is all that is important here. The transformation matrix is of the form:

$$
M = \begin{pmatrix} r_{11}&r_{12}&r_{13}&t_x\\ r_{21}&r_{22}&r_{23}&t_y \\ r_{31}&r_{32}&r_{33}&t_z \\ 0 & 0& 0& 1\end{pmatrix}
$$

Where the matrix $R \in SO(3)$ is a rotation matrix and $\vec t$ the translation vector.

## The algorithm

To construct the helix, we will assume that the first point B lies at the origin. 
The first segment goes in the x direction, so C is at $(L,0,0,1)$.
Further, the joint faces of the first segment are normal to the $xy$-plane, meaning
$$N_B = (-\cos(\theta/2), \sin(\theta/2), 0, 1)\quad N_C = (\cos(\theta/2), \sin(\theta/2), 0, 1)$$
which are the normal vectors to the joint faces.

Note that this means that the axis of the helix is not the $z$-axis, and it is in fact not even parallel to it unless the twist angle $\tau$ is 0.

To determine the third point A we do the following:

1. translate C to the origin ($\vec t = (-L, 0, 0), R = \mathbf 1$)
2. rotate by $-\theta/2$ along the $z$-axis. Now $\vec N_B$ is aligned with the $x$-axis.
3. rotate by $\tau$ around the $x$-axis.
4. rotate by $-\theta/2$ around the $z$-axis again. Now the point D would lie where C was before, and B has moved to A.

Doing this transformation to B will give you A, and doing the inverse to C will give you D.