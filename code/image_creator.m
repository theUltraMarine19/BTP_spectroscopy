% for minerals
function img = image_creator(num_channels)
img = zeros(300, 300, num_channels);
mineral1 = [9.1297, 9.130600000000001, 9.148, 9.2026, 9.1893, 9.1996, 9.2026, 9.1706, 8.875599999999999, 7.986000000000001, 6.6249, 5.324999999999999, 4.2878, 3.4939, 2.9031000000000002, 2.4645, 2.1422, 1.9064999999999999, 1.7323999999999997, 1.6077, 1.5185, 1.4529, 1.4058, 1.3727, 1.3533, 1.3527, 1.375, 1.4158, 1.4405999999999999, 1.4806, 1.4819, 1.4387, 1.4391, 1.4783, 1.5449000000000002, 1.6882, 1.7368000000000001, 1.8125, 1.9677, 1.9067, 1.795, 1.77, 1.7924, 1.7719, 1.8017999999999998, 1.9097999999999997, 1.9611, 2.2359, 2.6407, 2.7256, 2.4608999999999996, 2.1893000000000002, 1.9518000000000002, 1.7878, 1.7106, 1.7047999999999999, 1.7306, 1.7748, 1.8218, 1.8304999999999998, 1.8047, 1.7428, 1.7157999999999998, 1.704, 1.6794000000000002, 1.6418, 1.5935000000000001, 1.5464, 1.5083000000000002, 1.478, 1.4581, 1.4491, 1.4499, 1.4684, 1.5089000000000001, 1.5359, 1.505, 1.4264999999999999, 1.3374000000000001, 1.2545, 1.1785, 1.1197, 1.0779, 1.0493999999999999, 1.0282, 1.0123, 1.0004, 0.99, 0.9773000000000001, 0.9575, 0.9314, 0.9069, 0.8894, 0.8742, 0.8618, 0.8513000000000001, 0.8436999999999999, 0.8372999999999999, 0.8316999999999999, 0.8270000000000001, 0.8217000000000001, 0.8186, 0.8158000000000001, 0.8122, 0.8094, 0.8069, 0.8041, 0.8023000000000001, 0.7990999999999999, 0.7977000000000001, 0.7949, 0.7931999999999999, 0.7909999999999999, 0.7877000000000001, 0.7868999999999999, 0.7844, 0.7814000000000001, 0.7809999999999999, 0.7788999999999999, 0.7764, 0.7747999999999999, 0.7734000000000001, 0.7719999999999999, 0.7695000000000001, 0.768, 0.7661, 0.7646999999999999, 0.7629, 0.7602, 0.7587999999999999, 0.7571000000000001, 0.7543, 0.7519, 0.7508, 0.748, 0.7454, 0.7433000000000001, 0.7408, 0.7387, 0.7356, 0.733, 0.7307, 0.7281, 0.7241, 0.722, 0.7187, 0.7159, 0.712, 0.7099, 0.707, 0.7036, 0.7023, 0.6993, 0.6971999999999999, 0.6937, 0.6921999999999999, 0.6913, 0.6877, 0.6847, 0.6805, 0.6798, 0.6752, 0.6716, 0.6674, 0.6646000000000001, 0.6624, 0.6584, 0.6554, 0.6524, 0.6509, 0.6473, 0.6443000000000001, 0.6415, 0.6405000000000001, 0.6369, 0.6335, 0.6323, 0.6315999999999999, 0.6285000000000001, 0.6266, 0.6264, 0.6227, 0.6248, 0.6232000000000001, 0.6212, 0.6199, 0.6191, 0.6213, 0.6185, 0.6155999999999999, 0.6135999999999999, 0.6108, 0.6069];
mineral2 = [1.2494, 1.2111, 1.1783000000000001, 1.1531, 1.147, 1.1604, 1.3737, 1.0909, 0.978, 0.9146000000000001, 0.8733, 0.8407, 0.8142999999999999, 0.7917000000000001, 0.774, 0.761, 0.7533, 0.7524, 0.7589, 0.7742, 0.7990999999999999, 0.8321, 0.867, 0.8895, 0.8875, 0.8594, 0.8137, 0.7668999999999999, 0.7319, 0.7087, 0.6936, 0.6838, 0.6774, 0.6728, 0.6697, 0.6673, 0.6648000000000001, 0.6632, 0.6608, 0.6587000000000001, 0.656, 0.6537999999999999, 0.6513, 0.6488999999999999, 0.6469, 0.6448, 0.643, 0.642, 0.6412, 0.6402, 0.6396999999999999, 0.6396, 0.6397999999999999, 0.6397999999999999, 0.6401, 0.6406000000000001, 0.6419, 0.6426000000000001, 0.6437, 0.6443, 0.6454000000000001, 0.6469, 0.6483, 0.6493, 0.6497999999999999, 0.6482, 0.6447, 0.64, 0.6355999999999999, 0.6342, 0.6355, 0.6388, 0.6398, 0.6408, 0.6415, 0.6406000000000001, 0.6363, 0.6298, 0.6202, 0.6093, 0.5982999999999999, 0.5964, 0.5995999999999999, 0.6069, 0.6157, 0.6192, 0.6192, 0.618, 0.6143, 0.6112, 0.6081000000000001, 0.6053999999999999, 0.6045, 0.6038, 0.6055999999999999, 0.6077, 0.6109, 0.6166, 0.6243, 0.6335, 0.6453, 0.6595, 0.6757, 0.6948, 0.7165, 0.7404, 0.7678, 0.7979, 0.8316, 0.8679000000000001, 0.9075, 0.9507, 0.9977, 1.0473, 1.1011, 1.1571, 1.2173, 1.2809, 1.3476000000000001, 1.4192, 1.494, 1.573, 1.6557, 1.7436, 1.8353000000000002, 1.9317000000000002, 2.0326, 2.1375, 2.2449, 2.3563, 2.4721, 2.5933, 2.7184, 2.8472999999999997, 2.9816000000000003, 3.1181, 3.2576, 3.3998, 3.5443999999999996, 3.6910999999999996, 3.8395, 3.9908, 4.1409, 4.292, 4.4418999999999995, 4.5914, 4.7434, 4.896, 5.0437, 5.1894, 5.334099999999999, 5.4708000000000006, 5.6114999999999995, 5.7456, 5.8694, 5.995, 6.120000000000001, 6.2347, 6.339499999999999, 6.4267, 6.5091, 6.5894, 6.655200000000001, 6.7166, 6.7721, 6.8151, 6.8576, 6.890199999999999, 6.912599999999999, 6.9351, 6.9444, 6.9548, 6.9505, 6.9437, 6.9327, 6.9148, 6.8928, 6.8574, 6.811999999999999, 6.7654, 6.719799999999999, 6.6594999999999995, 6.5855999999999995, 6.5114, 6.4288, 6.3363, 6.2219999999999995, 6.0964, 5.9597, 5.8344000000000005, 5.7097, 5.5944, 5.4733];
mineral3 = [6.3032, 6.0941, 5.8871, 5.6854000000000005, 5.486400000000001, 5.3011, 5.1187000000000005, 4.9437, 4.7732, 4.6074, 4.4478, 4.2966, 4.1519, 4.0137, 3.8808999999999996, 3.7534, 3.6292, 3.5092, 3.3951000000000002, 3.2862, 3.1820000000000004, 3.0833, 2.9885, 2.8972, 2.8083, 2.7218999999999998, 2.638, 2.5573, 2.4813, 2.4104, 2.3444, 2.2819000000000003, 2.2218999999999998, 2.163, 2.107, 2.0522, 2.0002, 1.9518, 1.9068, 1.866, 1.8275000000000001, 1.7912, 1.7566, 1.7227999999999999, 1.6905000000000001, 1.6594, 1.6299, 1.6023, 1.5753, 1.5513, 1.5286, 1.5088999999999997, 1.4890999999999999, 1.4714, 1.4541, 1.4385, 1.4226999999999999, 1.4076, 1.3925, 1.3773, 1.3638, 1.3514, 1.3408, 1.3321, 1.3245, 1.3151, 1.3075, 1.3008, 1.2937, 1.2872000000000001, 1.2804, 1.2731999999999999, 1.2662, 1.2593, 1.2524, 1.2462, 1.2397, 1.2342, 1.2278, 1.2234, 1.2216, 1.2181, 1.2165, 1.2149999999999999, 1.2122, 1.2113, 1.2086000000000001, 1.2029, 1.1985, 1.1939, 1.1886999999999999, 1.1842, 1.1797, 1.1761, 1.1734, 1.1711, 1.1678, 1.165, 1.1615000000000002, 1.1583, 1.157, 1.1576, 1.1571, 1.155, 1.1524, 1.1502999999999999, 1.15, 1.1505, 1.1495, 1.1442, 1.1396, 1.1385, 1.1351, 1.1303999999999998, 1.1256, 1.1224, 1.1214, 1.1176, 1.1133, 1.113, 1.1137, 1.1109, 1.1061, 1.1084, 1.1088999999999998, 1.1041, 1.1047, 1.1058000000000001, 1.1026, 1.1009, 1.1028, 1.0996, 1.0967, 1.0984, 1.0948, 1.0911, 1.0937999999999999, 1.088, 1.0841, 1.0848, 1.0781999999999998, 1.0762, 1.0762, 1.0694, 1.0711, 1.0695000000000001, 1.067, 1.0679, 1.0651, 1.0678, 1.0665, 1.0665, 1.0701, 1.0652, 1.0685, 1.0661, 1.0629, 1.0657, 1.0608, 1.0636999999999999, 1.0564, 1.0609, 1.0561, 1.057, 1.0551, 1.0518999999999998, 1.0542, 1.0463, 1.0508, 1.0422, 1.0483, 1.0379, 1.0424, 1.0352999999999999, 1.0354, 1.0312, 1.0324, 1.0279, 1.0290000000000001, 1.0241, 1.0262, 1.0231, 1.0247, 1.0222000000000002, 1.0245, 1.0198, 1.0199, 1.0146, 1.0181, 1.0145, 1.0214, 1.0148, 1.0218];
mineral4 = [5.2242, 5.1766000000000005, 5.1290000000000004, 5.0714, 5.0221, 4.9681, 4.9126, 4.8574, 4.8211, 4.7863, 4.7646, 4.754, 4.7456, 4.7371, 4.7211, 4.7036, 4.6746, 4.6541, 4.6426, 4.6258, 4.6187000000000005, 4.6101, 4.5991, 4.5944, 4.590400000000001, 4.5943000000000005, 4.599, 4.6083, 4.6295, 4.6486, 4.6735, 4.6954, 4.719, 4.7357, 4.763, 4.7822, 4.804200000000001, 4.8261, 4.8458, 4.8687000000000005, 4.896599999999999, 4.9278, 4.9627, 4.9904, 5.0337, 5.0655, 5.1136, 5.1535, 5.2001, 5.246499999999999, 5.2937, 5.3261, 5.3846, 5.4469, 5.509, 5.5716, 5.6602, 5.698499999999999, 5.7223, 5.7524, 5.771, 5.8027, 5.8191, 5.8408999999999995, 5.8568, 5.8734, 5.8861, 5.9049000000000005, 5.92, 5.9383, 5.9593, 5.9755, 5.992, 6.0223, 6.0298, 6.0485999999999995, 6.0739, 6.085699999999999, 6.106999999999999, 6.109999999999999, 6.117, 6.1206, 6.1332, 6.1482, 6.144399999999999, 6.1472, 6.1459, 6.1443, 6.1567, 6.1457, 6.1447, 6.1467, 6.1533999999999995, 6.1538, 6.1541999999999994, 6.168699999999999, 6.190199999999999, 6.181900000000001, 6.1675, 6.160299999999999, 6.162600000000001, 6.1745, 6.1722, 6.165699999999999, 6.1492, 6.1404, 6.156699999999999, 6.141299999999999, 6.1504, 6.1581, 6.1505, 6.165900000000001, 6.1813, 6.1646, 6.1633, 6.1681, 6.1701999999999995, 6.1814, 6.1767, 6.174300000000001, 6.185499999999999, 6.1785, 6.1807, 6.194799999999999, 6.192500000000001, 6.1885, 6.177899999999999, 6.1882, 6.194799999999999, 6.2363, 6.246499999999999, 6.304600000000001, 6.3323, 6.389799999999999, 6.4046, 6.4873, 6.5977, 6.6029, 6.714600000000001, 6.718999999999999, 6.878500000000001, 6.9458, 6.9026, 6.986699999999999, 6.979900000000001, 6.963699999999999, 6.8339, 6.9373000000000005, 6.8304, 6.721200000000001, 6.6533999999999995, 6.6851, 6.4977, 6.4483999999999995, 6.331799999999999, 6.229100000000001, 6.187799999999999, 6.0736, 6.0321, 6.048500000000001, 6.0634, 6.0, 5.9732, 5.933, 5.8874, 5.8454, 5.914300000000001, 5.9166, 5.8277, 5.8222, 5.7682, 5.7818, 5.689, 5.5515, 5.5678, 5.6011, 5.5679, 5.547899999999999, 5.4952000000000005, 5.4779, 5.4107, 5.3424000000000005, 5.339700000000001, 5.2133, 5.2296000000000005, 5.3198, 5.5014, 5.9921, 7.0427, 36.351600000000005, 36.351600000000005, 36.351600000000005, 36.351600000000005];
mineral5 = [1.2285, 1.1817, 1.1378, 1.0996, 1.0638999999999998, 1.0326, 1.0043, 0.9783, 0.9556, 0.9346000000000001, 0.9159000000000002, 0.899, 0.8845000000000001, 0.8721000000000001, 0.8606, 0.8517, 0.8435999999999999, 0.8353999999999999, 0.8306, 0.8279000000000001, 0.8265, 0.8258, 0.827, 0.8285, 0.8313, 0.8343, 0.8375000000000001, 0.8423999999999999, 0.8478999999999999, 0.8526, 0.8574999999999999, 0.8627, 0.8683000000000001, 0.8744000000000001, 0.8802000000000001, 0.8848999999999999, 0.8917000000000002, 0.8969, 0.9019999999999999, 0.9062, 0.9103000000000001, 0.9136, 0.9187000000000001, 0.9242999999999999, 0.9295, 0.937, 0.9444, 0.9541, 0.9654, 0.9774, 0.9908999999999999, 1.0049000000000001, 1.022, 1.0392000000000001, 1.0594000000000001, 1.0801, 1.1014, 1.1239, 1.1487, 1.1760000000000002, 1.2045, 1.2333, 1.2634999999999998, 1.2949000000000002, 1.3279, 1.3624, 1.3977, 1.4346999999999999, 1.4721, 1.5114, 1.5499, 1.5903999999999998, 1.6306999999999998, 1.6730999999999998, 1.7155, 1.7585000000000002, 1.8034999999999999, 1.8465, 1.8908, 1.9369, 1.9808, 2.0278, 2.0738000000000003, 2.1212, 2.1668000000000003, 2.2154, 2.2633, 2.3099, 2.3587, 2.4067, 2.4554, 2.5008, 2.5502, 2.5968, 2.6448, 2.6913, 2.7389, 2.7861000000000002, 2.8325, 2.8791, 2.9255, 2.975, 3.0204, 3.0663, 3.1133999999999995, 3.1603000000000003, 3.2045999999999997, 3.2511, 3.2960000000000003, 3.3416000000000006, 3.3847000000000005, 3.4291, 3.4737, 3.5195000000000003, 3.5620000000000003, 3.6033000000000004, 3.6493, 3.6936999999999998, 3.738, 3.7814, 3.8253000000000004, 3.8712999999999997, 3.918, 3.9625000000000004, 4.0074, 4.053599999999999, 4.0988, 4.1463, 4.199999999999999, 4.245, 4.2923, 4.3374, 4.3832, 4.4261, 4.4776, 4.5137, 4.563599999999999, 4.609999999999999, 4.6543, 4.7044, 4.756399999999999, 4.8017, 4.8432, 4.8952, 4.9416, 4.9898, 5.0348999999999995, 5.0889, 5.1394, 5.180400000000001, 5.23, 5.285, 5.3299, 5.3688, 5.416399999999999, 5.4829, 5.5111, 5.5548, 5.614800000000001, 5.6437, 5.6888000000000005, 5.7579, 5.773, 5.808, 5.8408, 5.8891, 5.9331, 5.9459, 5.9893, 6.0319, 6.059000000000001, 6.077999999999999, 6.117500000000001, 6.1338, 6.1689, 6.182, 6.1973, 6.229800000000001, 6.2347, 6.2441, 6.274100000000001, 6.255199999999999, 6.3079, 6.2706, 6.311999999999999, 6.2867999999999995, 6.2705, 6.2607, 6.2447, 6.2432, 6.2509, 6.2132000000000005, 6.2394];

x_cens = randi([0 200], 1, 4);
y_cens = randi([0 200], 1, 4);
widths = randi([50 100], 1, 4);
lengths = randi([50 100], 1, 4);
for i=1:num_channels
    img(:,:,i) = mineral2(i);
    img(x_cens(1):x_cens(1)+lengths(1), y_cens(1):y_cens(1)+widths(1),i) = mineral1(i);
    img(x_cens(2):x_cens(2)+lengths(2), y_cens(2):y_cens(2)+widths(2),i) = mineral3(i);
    img(x_cens(3):x_cens(3)+lengths(3), y_cens(3):y_cens(3)+widths(3),i) = mineral4(i);
    img(x_cens(4):x_cens(4)+lengths(4), y_cens(4):y_cens(4)+widths(4),i) = mineral5(i);
end