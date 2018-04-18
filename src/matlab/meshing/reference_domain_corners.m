function C = reference_domain_corners(id)
switch id
    case ShapeSet.LinearL
        C = [
            0 0
            1 0
            0 1
            ];
    case 24
        C = [
            -1 -1
            +1 -1
            +1 +1
            -1 +1
            ];
end
end