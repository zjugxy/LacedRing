1.所有face初始化为dualnodetype

faceidx --> meshletidx --> dualtype

# tradeoff
\
To generate meshlet data, this library provides 
two algorithms - meshopt_buildMeshletsScan, which creates
the meshlet data using a vertex cache-optimized index buffer 
as a starting point by greedily aggregating consecutive triangles
until they go over the meshlet limits, and meshopt_buildMeshlets,
which doesn't depend on any other algorithms and tries to balance
topological efficiency (by maximizing vertex reuse inside meshlets)
with culling efficiency (by minimizing meshlet radius and triangle
direction divergence). meshopt_buildMeshlets is recommended in cases
when the resulting meshlet data 
will be used in cluster culling algorithms.

https://github.com/zeux/meshoptimizer#mesh-shading

可以研究一下不同种类meshlet对渲染和压缩的影响。