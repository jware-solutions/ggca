mod common_adjustment;

use crate::common_adjustment::test_adjustment_bh_or_by;
use ggca::adjustment::AdjustmentMethod;

#[test]
fn test_benjamini_yekutieli_1() {
    let x = vec![
        0.492476438,
        0.581083727,
        0.210531377,
        0.789986095,
        0.589752191,
        0.264190242,
        0.494516471,
        0.606625965,
        0.475790917,
        0.645544453,
        0.297524352,
        0.011062263,
        0.010853499,
        0.117684559,
        0.457888339,
        0.117206734,
        0.138873356,
        0.109990668,
        0.512104486,
        0.220523547,
        0.172112147,
        0.650634463,
        0.739423803,
        0.217546492,
        0.931440274,
        0.333091045,
        0.897488627,
        0.462819061,
        0.380017506,
        0.139116981,
        0.098543480,
        0.852069447,
        0.334188761,
        0.589486521,
        0.089390080,
        0.038774872,
        0.798862100,
        0.864838914,
        0.167394274,
        0.149253015,
        0.949378401,
        0.610268435,
        0.744387908,
        0.597339348,
        0.692591534,
        0.490572178,
        0.587196384,
        0.005123415,
        0.438803020,
        0.200521362,
        0.216752211,
        0.812404521,
        0.774398944,
        0.301477611,
        0.290656052,
        0.369112384,
        0.941408090,
        0.623077547,
        0.289655600,
        0.010498420,
        0.513816550,
        0.935371934,
        0.388870449,
        0.579237843,
        0.607088043,
        0.062218425,
        0.502608138,
        0.984300127,
        0.700455081,
        0.710313653,
        0.305986631,
        0.429015495,
        0.336344790,
        0.035646178,
        0.309525103,
        0.449058992,
        0.262996294,
        0.700174163,
        0.689939301,
        0.683347777,
        0.875076623,
        0.664084711,
        0.865176424,
        0.026736447,
        0.140241586,
        0.024521381,
        0.954727222,
        0.186416370,
        0.885015741,
        0.893631153,
        0.739436367,
        0.596530028,
        0.516476509,
        0.378430589,
        0.848003399,
        0.102478122,
        0.226252391,
        0.462288224,
        0.639915368,
        0.072198368,
    ];

    let expected = vec![
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    ];

    test_adjustment_bh_or_by(x, expected, AdjustmentMethod::BenjaminiYekutieli);
}

#[test]
fn test_benjamini_yekutieli_2() {
    let x = vec![
        5.39198144708e-77,
        1.18917839281e-160,
        9.48402847878e-209,
        3.82935778527e-82,
        2.67977253966e-114,
    ];

    let expected = vec![
        1.23116910e-076,
        6.78822666e-160,
        1.08275992e-207,
        1.09296253e-081,
        1.01980233e-113,
    ];

    test_adjustment_bh_or_by(x, expected, AdjustmentMethod::BenjaminiYekutieli);
}

#[test]
fn test_benjamini_yekutieli_3() {
    let x = vec![
        5.455469E-11,
        8.966458E-11,
        4.444592E-10,
        8.539654E-10,
        1.336221E-09,
        2.087884E-09,
        2.216075E-09,
        2.331999E-09,
        3.207092E-09,
        4.534499E-09,
        4.619085E-09,
        5.827716E-09,
    ];

    let expected = vec![
        1.66948849e-09,
        1.66948849e-09,
        5.51700214e-09,
        7.95010364e-09,
        9.95178066e-09,
        1.08550263e-08,
        1.08550263e-08,
        1.08550263e-08,
        1.32697095e-08,
        1.56370842e-08,
        1.56370842e-08,
        1.80846305e-08,
    ];

    test_adjustment_bh_or_by(x, expected, AdjustmentMethod::BenjaminiYekutieli);
}