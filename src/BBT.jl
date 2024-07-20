####################################################
# BBT - BoundingBoxTree
#
# A tree that stores points in a hierarchical axis aligned bounding boxes.
####################################################

module  BBT

IRange = UnitRange{Int64};
mutable struct BBTNode
    bb::BBox{D,T}
    r::IRnage
    left::Union{Nothing, BBTNode}
    right::Union{Nothing, BBTNode}
end

mutable struct BoundingBoxTree{D,T}
    PS::Polygon{D,T}
    root::BBTNode 
end

function  node_init( PS::Polygon{D,T}, range::IRange )
    bb = BBox_init( PS, range );
    node = BBTNode( bb, range, Nothing, Nothing );
    return  node;
end


function  node_split( node::BBTNode, PS::Polygon{D,T} )
    ( length( node.r ) <= 1 )  return;

    
    
    return  node;
end
