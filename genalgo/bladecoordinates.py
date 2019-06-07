import numpy as np
import matplotlib.pyplot as plt
blade2=np.array([[0.1,0.200079,0],
[0.098418996,0.200349706,0],
[0.096880033,0.200605949,0],
[0.095374971,0.200846456,0],
[0.09381446,0.201086868,0],
[0.091652794,0.201406068,0],
[0.090068341,0.201630028,0],
[0.088229058,0.20187954,0],
[0.086827398,0.202062068,0],
[0.08598648,0.202168376,0],
[0.083675243,0.202448001,0],
[0.082170866,0.20261968,0],
[0.081371922,0.202707356,0],
[0.080660503,0.202783311,0],
[0.07905283,0.202947333,0],
[0.077501637,0.203094916,0],
[0.076129225,0.203216007,0],
[0.074752211,0.203327895,0],
[0.073333744,0.20343226,0],
[0.071875391,0.20352676,0],
[0.070498554,0.203602833,0],
[0.069167449,0.203663106,0],
[0.067431717,0.203719759,0],
[0.064105822,0.203746562,0],
[0.06220189,0.203703465,0],
[0.060291794,0.203606274,0],
[0.058523854,0.203458113,0],
[0.057585053,0.203351714,0],
[0.05680216,0.203245692,0],
[0.055752435,0.203074053,0],
[0.054838609,0.202891853,0],
[0.054361066,0.202782362,0],
[0.053773599,0.202631475,0],
[0.053150983,0.202448001,0],
[0.052532709,0.202233895,0],
[0.052105547,0.202062068,0],
[0.051694359,0.201872795,0],
[0.051243226,0.201630028,0],
[0.050840391,0.201364189,0],
[0.050526563,0.201086868,0],
[0.050168866,0.200605949,0],
[0.050058204,0.200349706,0],
[0.050003115,0.200079,0],
[0.050061086,0.199641499,0],
[0.050223237,0.19929986,0],
[0.050526563,0.198913132,0],
[0.050840391,0.198635811,0],
[0.051694359,0.198127205,0],
[0.051243226,0.198369972,0],
[0.052532709,0.197766105,0],
[0.053150983,0.197551999,0],
[0.053773599,0.197368525,0],
[0.054361066,0.197217638,0],
[0.054838609,0.197108147,0],
[0.055752435,0.196925947,0],
[0.056332189,0.196826621,0],
[0.05680216,0.196754308,0],
[0.057585053,0.196648286,0],
[0.058523854,0.196541887,0],
[0.060291794,0.196393726,0],
[0.06220189,0.196296535,0],
[0.064105822,0.196253438,0],
[0.065614912,0.196251139,0],
[0.067431717,0.196280241,0],
[0.069167449,0.196336894,0],
[0.070498554,0.196397167,0],
[0.071875391,0.19647324,0],
[0.073333744,0.19656774,0],
[0.077501637,0.196905084,0],
[0.07905283,0.197052667,0],
[0.080660503,0.197216689,0],
[0.081371922,0.197292644,0],
[0.082170866,0.19738032,0],
[0.08598648,0.197831624,0],
[0.083675243,0.197551999,0],
[0.08598648,0.197831624,0],
[0.086827398,0.197937932,0],
[0.088229058,0.19812046,0],
[0.088229058,0.19812046,0],
[0.090068341,0.198369972,0],
[0.091652794,0.198593932,0],
[0.09381446,0.198913132,0],
[0.095374971,0.199153544,0],
[0.096880033,0.199394051,0],
[0.098418996,0.199650294,0],
[0.1,0.200079,0]])




blade1=np.array([[0.036127975,0.20971822,0],
[0.036167406,0.209908331,0],
[0.034590653,0.209756988,0],
[0.033468245,0.209645186,0],
[0.032178942,0.209509938,0],
[0.030615263,0.209337112,0],
[0.029014409,0.209151435,0],
[0.027347488,0.208948919,0],
[0.025631702,0.208730815,0],
[0.023874121,0.208497268,0],
[0.021667409,0.208189244,0],
[0.019498395,0.207869979,0],
[0.017716674,0.207594684,0],
[0.016187801,0.20734839,0],
[0.01456403,0.207075913,0],
[0.013193908,0.206836498,0],
[0.011432248,0.206514688,0],
[0.009767676,0.206194792,0],
[0.008374207,0.205913655,0],
[0.005683873,0.205331042,0],
[0.006527214,0.205519769,0],
[0.004844565,0.20513714,0],
[0.003077,0.204706773,0],
[0.001364009,0.204257557,0],
[-0.00030658,0.203784041,0],
[-0.002163435,0.203207584,0],
[-0.003670457,0.202692556,0],
[-0.004610316,0.202345704,0],
[-0.005037683,0.202180366,0],
[-0.005739721,0.201897181,0],
[-0.006408675,0.201612388,0],
[-0.007000892,0.201346013,0],
[-0.007571615,0.201074628,0],
[-0.006122774,0.201736044,0],
[-0.008026844,0.200846456,0],
[-0.008483054,0.200605949,0],
[-0.009063441,0.200279412,0],
[-0.009397488,0.200079,0],
[-0.009724868,0.199871495,0],
[-0.010065084,0.199641499,0],
[-0.010379474,0.199413117,0],
[-0.010448592,0.199360527,0],
[-0.010594411,0.199246348,0],
[-0.010526793,0.19929986,0],
[-0.010682875,0.19917476,0],
[-0.010802344,0.199075006,0],
[-0.010897384,0.198992889,0],
[-0.010986848,0.198913132,0],
[-0.011071864,0.198834932,0],
[-0.011138796,0.198771574,0],
[-0.011220114,0.198692289,0],
[-0.011317251,0.198593932,0],
[-0.011438872,0.198463604,0],
[-0.011520314,0.198369972,0],
[-0.011645665,0.198211366,0],
[-0.011751312,0.198058039,0],
[-0.011843939,0.197901484,0],
[-0.011948389,0.197685768,0],
[-0.011912529,0.197766105,0],
[-0.012000152,0.197551999,0],
[-0.012026936,0.1974699,0],
[-0.012054612,0.197368525,0],
[-0.012071277,0.197292644,0],
[-0.012096259,0.197108147,0],
[-0.012096259,0.196920093,0],
[-0.012099511,0.1969883,0],
[-0.012085152,0.196826621,0],
[-0.012054612,0.196694492,0],
[-0.012023014,0.196605364,0],
[-0.011817641,0.196268748,0],
[-0.011520314,0.195982975,0],
[-0.011994893,0.196541887,0],
[-0.011959099,0.19647324,0],
[-0.011910446,0.196393726,0],
[-0.01170119,0.196141628,0],
[-0.011624205,0.196069468,0],
[-0.011401829,0.195896349,0],
[-0.011138796,0.195738926,0],
[-0.011273337,0.195814082,0],
[-0.010986848,0.195665826,0],
[-0.010802344,0.19559111,0],
[-0.010682875,0.19554943,0],
[-0.010526793,0.195501195,0],
[-0.010379474,0.19546059,0],
[-0.010205376,0.195417344,0],
[-0.010065084,0.195385906,0],
[-0.009909973,0.195354433,0],
[-0.00974637,0.195324717,0],
[-0.009582132,0.195298223,0],
[-0.009397488,0.195272131,0],
[-0.009240227,0.195252767,0],
[-0.00910465,0.195238037,0],
[-0.009001581,0.195227987,0],
[-0.008849148,0.195214855,0],
[-0.008691602,0.195203361,0],
[-0.008564088,0.19519553,0],
[-0.008400002,0.195187291,0],
[-0.008223946,0.195180621,0],
[-0.008026844,0.195175635,0],
[-0.007571615,0.195173005,0],
[-0.007842846,0.195173163,0],
[-0.007257558,0.195177618,0],
[-0.007000892,0.195184966,0],
[-0.006807554,0.195192496,0],
[-0.006632197,0.19520073,0],
[-0.006408675,0.195213061,0],
[-0.006251679,0.195222898,0],
[-0.006122774,0.195231683,0],
[-0.005989557,0.195241419,0],
[-0.005810353,0.195255548,0],
[-0.005676272,0.195266876,0],
[-0.005409527,0.195291287,0],
[-0.005037683,0.195329312,0],
[-0.004549277,0.195385906,0],
[-0.003980275,0.19546059,0],
[-0.003670457,0.195504889,0],
[-0.003020217,0.195605642,0],
[-0.002477769,0.195697314,0],
[-0.001937655,0.195795011,0],
[-0.001565879,0.195865751,0],
[-0.001187306,0.19594054,0],
[-0.000735364,0.196033329,0],
[-0.00030658,0.196124775,0],
[0.000188707,0.196234384,0],
[0.000747092,0.196362861,0],
[0.001364009,0.196510547,0],
[0.001747954,0.196605364,0],
[0.002181115,0.196714899,0],
[0.002612655,0.196826621,0],
[0.003077,0.196949607,0],
[0.003765258,0.197136901,0],
[0.004263198,0.197275966,0],
[0.004506739,0.197345039,0],
[0.004822186,0.197435518,0],
[0.005089566,0.197513089,0],
[0.005683873,0.19768832,0],
[0.006527214,0.197943392,0],
[0.008374207,0.198526295,0],
[0.009767676,0.198985977,0],
[0.011432248,0.199555304,0],
[0.01448996,0.200651066,0],
[0.013193908,0.200179334,0],
[0.016187801,0.201283982,0],
[0.017716674,0.201867338,0],
[0.019498395,0.202561925,0],
[0.021667409,0.203427569,0],
[0.023874121,0.204329356,0],
[0.025631702,0.205062279,0],
[0.026713429,0.205519769,0],
[0.027794219,0.205981702,0],
[0.029014409,0.206509083,0],
[0.030615263,0.207210611,0],
[0.032178942,0.207906636,0],
[0.033752568,0.208619225,0],
[0.034590653,0.209004538,0],
[0.035084788,0.209233271,0],
[0.035495522,0.209423983,0],
[0.035814419,0.209572289,0],
[0.036127975,0.20971822,0]])





# plt.plot(blade2.transpose()[0],blade2.transpose()[1],'.')
# plt.plot(blade1.transpose()[0],blade1.transpose()[1],'.')
# plt.axes().set_aspect('equal')
# plt.show()
# [0.050133367,0.199463714,0],
# [0.050359597,0.199104668,0],
# [0.050692087,0.198756158,0],
# [0.051020536,0.198508325,0],
# [0.05145099,0.198252931,0],
# [0.052105547,0.197937932,0],
# [0.052230575,0.197885291,0],
# [0.051920551,0.19801983,0],
# [0.051770774,0.198089978,0],
# [0.050008164,0.199871495,0],
# [0.050032876,0.199739053,0],
# [0.05000004,0.199991076,0],
# [0.050028816,0.200243952,0],
# [0.050092656,0.200444301,0],
# [0.050133367,0.200536286,0],
# [0.050238753,0.200724904,0],
# [0.050193994,0.200651066,0],
# [0.050277114,0.200782926,0],
# [0.050322389,0.200846456,0],
# [0.050359597,0.200895332,0],
# [0.050418898,0.200968048,0],
# [0.050478345,0.201035537,0],
# [0.05001575,0.200179334,0],
# [0.050306885,0.19917476,0],
# [0.050418898,0.199031952,0],
# [0.050605403,0.198834932,0],
# [0.050948871,0.198556991,0],
# [0.051566793,0.198191686,0],
# [0.052829517,0.197658758,0],
# [0.052363249,0.197831624,0],
# [0.053368438,0.197484729,0],
# [0.055081097,0.197108147,0],
# [0.050092816,0.199555304,0],
# [0.050168866,0.199394051,0],
# [0.050605403,0.201165068,0],
# [0.050692087,0.201243842,0],
# [0.050948871,0.201443009,0],
# [0.051020536,0.201491675,0],
# [0.051146672,0.201572067,0],
# [0.05077746,0.201315105,0],
# [0.05135129,0.20169211,0],
# [0.051495846,0.201771112,0],
# [0.050007214,0.200120707,0],
# [0.050043737,0.200302002,0],
# [0.050071546,0.200388869,0],
# [0.050035197,0.200270211,0],
# [0.050021964,0.200212399,0],
# [0.050001382,0.199947521,0],
# [0.050003042,0.199921941,0],
# [0.050013716,0.199832833,0],
# [0.050021964,0.199787601,0],
# [0.050025992,0.199768556,0],
# [0.050043651,0.199698302,0],
# [0.050051665,0.199671057,0],
# [0.050000741,0.200038385,0],
# [0.050000225,0.200021103,0],
# [0.05001137,0.200151984,0],
# [0.050050693,0.200325752,0],
# [0.050149742,0.200569402,0],
# [0.050119277,0.200506201,0],
# [0.050102691,0.20046852,0],
# [0.050080224,0.200412489,0],
# [0.050064912,0.200369875,0],
# [0.050282088,0.199209852,0],
# [0.050257406,0.199246348,0],
# [0.050202941,0.199333559,0],
# [0.050187377,0.199360527,0],
# [0.050149742,0.199430598,0],
# [0.050108856,0.199517156,0],
# [0.050075983,0.199598894,0],
# [0.050333536,0.19913861,0],
# [0.050478345,0.198964463,0],
# [0.050395951,0.199075006,0],
# [0.050439807,0.199007645,0],
# [0.050016771,0.19981485,0],
# [0.050004719,0.200097403,0],
# [0.05000043,0.199970801,0],
# [0.050238753,0.199275096,0],
# [0.051146672,0.198427933,0],
# [0.0510713,0.198475239,0],
# [0.051324107,0.198323253,0],
# [0.051511778,0.198220447,0],
# [0.051393111,0.198284574,0],
# [0.051286575,0.198344742,0],
# [0.050768278,0.198692289,0],
# [0.050650021,0.198793509,0],
# [0.050896869,0.198593932,0],
# [0.050563077,0.19887608,0],
# [0.051920551,0.20198017,0],
# [0.051814065,0.201930671,0],
# [0.052033702,0.202030852,0],
# [0.052197238,0.202100875,0],
# [0.052283281,0.20213629,0],
# [0.052335961,0.202157514,0],
# [0.052391774,0.202179635,0],
# [0.052455163,0.202204317,0],
# [0.052598408,0.202258442,0],
# [0.051584028,0.201817203,0],
# [0.050005368,0.199896044,0],
# [0.084995461,0.202290552,0],
# [0.084995461,0.197709448,0],
# [0.1,0.199921,0],
# [0.09922372,0.200212399,0],
# [0.09922372,0.199787601,0],
# [0.099712206,0.199871495,0],
# [0.099712206,0.200128505,0],
# [0.1,0.200004183,0],
# [-0.011751312,0.196193255,0],
# [-0.011582416,0.196033329,0],
# [-0.010239121,0.199517156,0],
# [-0.009909973,0.199748367,0],
# [-0.009582132,0.199963471,0],
# [-0.007413809,0.195174701,0],
# [-0.005183624,0.195313848,0]])