/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
const
double VCT_L191[] = {
   /*     0  */  0.0000000E+00,
   /*     1  */  0.1989182E+01,
   /*     2  */  0.3978941E+01,
   /*     3  */  0.5969277E+01,
   /*     4  */  0.7960190E+01,
   /*     5  */  0.9951681E+01,
   /*     6  */  0.1194375E+02,
   /*     7  */  0.1393640E+02,
   /*     8  */  0.1592962E+02,
   /*     9  */  0.1792342E+02,
   /*    10  */  0.1991780E+02,
   /*    11  */  0.2169803E+02,
   /*    12  */  0.2346212E+02,
   /*    13  */  0.2522733E+02,
   /*    14  */  0.2700198E+02,
   /*    15  */  0.2879443E+02,
   /*    16  */  0.3061303E+02,
   /*    17  */  0.3246615E+02,
   /*    18  */  0.3436218E+02,
   /*    19  */  0.3630951E+02,
   /*    20  */  0.3831656E+02,
   /*    21  */  0.4040030E+02,
   /*    22  */  0.4255329E+02,
   /*    23  */  0.4477730E+02,
   /*    24  */  0.4707262E+02,
   /*    25  */  0.4945101E+02,
   /*    26  */  0.5192537E+02,
   /*    27  */  0.5450990E+02,
   /*    28  */  0.5722043E+02,
   /*    29  */  0.6006341E+02,
   /*    30  */  0.6304558E+02,
   /*    31  */  0.6617471E+02,
   /*    32  */  0.6945839E+02,
   /*    33  */  0.7290482E+02,
   /*    34  */  0.7652187E+02,
   /*    35  */  0.8031786E+02,
   /*    36  */  0.8430153E+02,
   /*    37  */  0.8848209E+02,
   /*    38  */  0.9286930E+02,
   /*    39  */  0.9747336E+02,
   /*    40  */  0.1023051E+03,
   /*    41  */  0.1073764E+03,
   /*    42  */  0.1126994E+03,
   /*    43  */  0.1182874E+03,
   /*    44  */  0.1241543E+03,
   /*    45  */  0.1303102E+03,
   /*    46  */  0.1367707E+03,
   /*    47  */  0.1435521E+03,
   /*    48  */  0.1506724E+03,
   /*    49  */  0.1581446E+03,
   /*    50  */  0.1659855E+03,
   /*    51  */  0.1742163E+03,
   /*    52  */  0.1828540E+03,
   /*    53  */  0.1919227E+03,
   /*    54  */  0.2014427E+03,
   /*    55  */  0.2114318E+03,
   /*    56  */  0.2219189E+03,
   /*    57  */  0.2329250E+03,
   /*    58  */  0.2444784E+03,
   /*    59  */  0.2566022E+03,
   /*    60  */  0.2693276E+03,
   /*    61  */  0.2826821E+03,
   /*    62  */  0.2967018E+03,
   /*    63  */  0.3114205E+03,
   /*    64  */  0.3268631E+03,
   /*    65  */  0.3430727E+03,
   /*    66  */  0.3600884E+03,
   /*    67  */  0.3779500E+03,
   /*    68  */  0.3966980E+03,
   /*    69  */  0.4163636E+03,
   /*    70  */  0.4370087E+03,
   /*    71  */  0.4586880E+03,
   /*    72  */  0.4814459E+03,
   /*    73  */  0.5053349E+03,
   /*    74  */  0.5304061E+03,
   /*    75  */  0.5567164E+03,
   /*    76  */  0.5843365E+03,
   /*    77  */  0.6133445E+03,
   /*    78  */  0.6437992E+03,
   /*    79  */  0.6757712E+03,
   /*    80  */  0.7093336E+03,
   /*    81  */  0.7445564E+03,
   /*    82  */  0.7815272E+03,
   /*    83  */  0.8203336E+03,
   /*    84  */  0.8610668E+03,
   /*    85  */  0.9038228E+03,
   /*    86  */  0.9487021E+03,
   /*    87  */  0.9958097E+03,
   /*    88  */  0.1045256E+04,
   /*    89  */  0.1097158E+04,
   /*    90  */  0.1151637E+04,
   /*    91  */  0.1208822E+04,
   /*    92  */  0.1268845E+04,
   /*    93  */  0.1331849E+04,
   /*    94  */  0.1397981E+04,
   /*    95  */  0.1467398E+04,
   /*    96  */  0.1540261E+04,
   /*    97  */  0.1616741E+04,
   /*    98  */  0.1697019E+04,
   /*    99  */  0.1781285E+04,
   /*   100  */  0.1869734E+04,
   /*   101  */  0.1962575E+04,
   /*   102  */  0.2060026E+04,
   /*   103  */  0.2162316E+04,
   /*   104  */  0.2269685E+04,
   /*   105  */  0.2382386E+04,
   /*   106  */  0.2500682E+04,
   /*   107  */  0.2624852E+04,
   /*   108  */  0.2755189E+04,
   /*   109  */  0.2891997E+04,
   /*   110  */  0.3035597E+04,
   /*   111  */  0.3186328E+04,
   /*   112  */  0.3344622E+04,
   /*   113  */  0.3510779E+04,
   /*   114  */  0.3685152E+04,
   /*   115  */  0.3868165E+04,
   /*   116  */  0.4060295E+04,
   /*   117  */  0.4261875E+04,
   /*   118  */  0.4473389E+04,
   /*   119  */  0.4695541E+04,
   /*   120  */  0.4928601E+04,
   /*   121  */  0.5173514E+04,
   /*   122  */  0.5430807E+04,
   /*   123  */  0.5700842E+04,
   /*   124  */  0.5984263E+04,
   /*   125  */  0.6281830E+04,
   /*   126  */  0.6593892E+04,
   /*   127  */  0.6921389E+04,
   /*   128  */  0.7265031E+04,
   /*   129  */  0.7624853E+04,
   /*   130  */  0.8001610E+04,
   /*   131  */  0.8396169E+04,
   /*   132  */  0.8808917E+04,
   /*   133  */  0.9240285E+04,
   /*   134  */  0.9689937E+04,
   /*   135  */  0.1015810E+05,
   /*   136  */  0.1064388E+05,
   /*   137  */  0.1114578E+05,
   /*   138  */  0.1166126E+05,
   /*   139  */  0.1218736E+05,
   /*   140  */  0.1272024E+05,
   /*   141  */  0.1325880E+05,
   /*   142  */  0.1380106E+05,
   /*   143  */  0.1434555E+05,
   /*   144  */  0.1488822E+05,
   /*   145  */  0.1542401E+05,
   /*   146  */  0.1594836E+05,
   /*   147  */  0.1645946E+05,
   /*   148  */  0.1695481E+05,
   /*   149  */  0.1743284E+05,
   /*   150  */  0.1788861E+05,
   /*   151  */  0.1831755E+05,
   /*   152  */  0.1871637E+05,
   /*   153  */  0.1908261E+05,
   /*   154  */  0.1941345E+05,
   /*   155  */  0.1970464E+05,
   /*   156  */  0.1995197E+05,
   /*   157  */  0.2015255E+05,
   /*   158  */  0.2030349E+05,
   /*   159  */  0.2040116E+05,
   /*   160  */  0.2044152E+05,
   /*   161  */  0.2042110E+05,
   /*   162  */  0.2033644E+05,
   /*   163  */  0.2018368E+05,
   /*   164  */  0.1995849E+05,
   /*   165  */  0.1965700E+05,
   /*   166  */  0.1927563E+05,
   /*   167  */  0.1881073E+05,
   /*   168  */  0.1825982E+05,
   /*   169  */  0.1761487E+05,
   /*   170  */  0.1687312E+05,
   /*   171  */  0.1603171E+05,
   /*   172  */  0.1508720E+05,
   /*   173  */  0.1403435E+05,
   /*   174  */  0.1286901E+05,
   /*   175  */  0.1159188E+05,
   /*   176  */  0.1019584E+05,
   /*   177  */  0.8788987E+04,
   /*   178  */  0.7438805E+04,
   /*   179  */  0.6144316E+04,
   /*   180  */  0.4941777E+04,
   /*   181  */  0.3850913E+04,
   /*   182  */  0.2887697E+04,
   /*   183  */  0.2063780E+04,
   /*   184  */  0.1385913E+04,
   /*   185  */  0.8553618E+03,
   /*   186  */  0.4673335E+03,
   /*   187  */  0.2103939E+03,
   /*   188  */  0.6588920E+02,
   /*   189  */  0.7367700E+01,
   /*   190  */  0.0000000E+00,
   /*   191  */  0.0000000E+00,

   /*     0  */  0.0000000E+00,
   /*     1  */  0.0000000E+00,
   /*     2  */  0.0000000E+00,
   /*     3  */  0.0000000E+00,
   /*     4  */  0.0000000E+00,
   /*     5  */  0.0000000E+00,
   /*     6  */  0.0000000E+00,
   /*     7  */  0.0000000E+00,
   /*     8  */  0.0000000E+00,
   /*     9  */  0.0000000E+00,
   /*    10  */  0.0000000E+00,
   /*    11  */  0.0000000E+00,
   /*    12  */  0.0000000E+00,
   /*    13  */  0.0000000E+00,
   /*    14  */  0.0000000E+00,
   /*    15  */  0.0000000E+00,
   /*    16  */  0.0000000E+00,
   /*    17  */  0.0000000E+00,
   /*    18  */  0.0000000E+00,
   /*    19  */  0.0000000E+00,
   /*    20  */  0.0000000E+00,
   /*    21  */  0.0000000E+00,
   /*    22  */  0.0000000E+00,
   /*    23  */  0.0000000E+00,
   /*    24  */  0.0000000E+00,
   /*    25  */  0.0000000E+00,
   /*    26  */  0.0000000E+00,
   /*    27  */  0.0000000E+00,
   /*    28  */  0.0000000E+00,
   /*    29  */  0.0000000E+00,
   /*    30  */  0.0000000E+00,
   /*    31  */  0.0000000E+00,
   /*    32  */  0.0000000E+00,
   /*    33  */  0.0000000E+00,
   /*    34  */  0.0000000E+00,
   /*    35  */  0.0000000E+00,
   /*    36  */  0.0000000E+00,
   /*    37  */  0.0000000E+00,
   /*    38  */  0.0000000E+00,
   /*    39  */  0.0000000E+00,
   /*    40  */  0.0000000E+00,
   /*    41  */  0.0000000E+00,
   /*    42  */  0.0000000E+00,
   /*    43  */  0.0000000E+00,
   /*    44  */  0.0000000E+00,
   /*    45  */  0.0000000E+00,
   /*    46  */  0.0000000E+00,
   /*    47  */  0.0000000E+00,
   /*    48  */  0.0000000E+00,
   /*    49  */  0.0000000E+00,
   /*    50  */  0.0000000E+00,
   /*    51  */  0.0000000E+00,
   /*    52  */  0.0000000E+00,
   /*    53  */  0.0000000E+00,
   /*    54  */  0.0000000E+00,
   /*    55  */  0.0000000E+00,
   /*    56  */  0.0000000E+00,
   /*    57  */  0.0000000E+00,
   /*    58  */  0.0000000E+00,
   /*    59  */  0.0000000E+00,
   /*    60  */  0.0000000E+00,
   /*    61  */  0.0000000E+00,
   /*    62  */  0.0000000E+00,
   /*    63  */  0.0000000E+00,
   /*    64  */  0.0000000E+00,
   /*    65  */  0.0000000E+00,
   /*    66  */  0.0000000E+00,
   /*    67  */  0.0000000E+00,
   /*    68  */  0.0000000E+00,
   /*    69  */  0.0000000E+00,
   /*    70  */  0.0000000E+00,
   /*    71  */  0.0000000E+00,
   /*    72  */  0.0000000E+00,
   /*    73  */  0.0000000E+00,
   /*    74  */  0.0000000E+00,
   /*    75  */  0.0000000E+00,
   /*    76  */  0.0000000E+00,
   /*    77  */  0.0000000E+00,
   /*    78  */  0.0000000E+00,
   /*    79  */  0.0000000E+00,
   /*    80  */  0.0000000E+00,
   /*    81  */  0.0000000E+00,
   /*    82  */  0.0000000E+00,
   /*    83  */  0.0000000E+00,
   /*    84  */  0.0000000E+00,
   /*    85  */  0.0000000E+00,
   /*    86  */  0.0000000E+00,
   /*    87  */  0.0000000E+00,
   /*    88  */  0.0000000E+00,
   /*    89  */  0.0000000E+00,
   /*    90  */  0.0000000E+00,
   /*    91  */  0.0000000E+00,
   /*    92  */  0.0000000E+00,
   /*    93  */  0.0000000E+00,
   /*    94  */  0.0000000E+00,
   /*    95  */  0.0000000E+00,
   /*    96  */  0.0000000E+00,
   /*    97  */  0.0000000E+00,
   /*    98  */  0.0000000E+00,
   /*    99  */  0.0000000E+00,
   /*   100  */  0.0000000E+00,
   /*   101  */  0.0000000E+00,
   /*   102  */  0.0000000E+00,
   /*   103  */  0.0000000E+00,
   /*   104  */  0.0000000E+00,
   /*   105  */  0.0000000E+00,
   /*   106  */  0.0000000E+00,
   /*   107  */  0.0000000E+00,
   /*   108  */  0.0000000E+00,
   /*   109  */  0.0000000E+00,
   /*   110  */  0.0000000E+00,
   /*   111  */  0.0000000E+00,
   /*   112  */  0.0000000E+00,
   /*   113  */  0.0000000E+00,
   /*   114  */  0.0000000E+00,
   /*   115  */  0.0000000E+00,
   /*   116  */  0.0000000E+00,
   /*   117  */  0.0000000E+00,
   /*   118  */  0.0000000E+00,
   /*   119  */  0.0000000E+00,
   /*   120  */  0.0000000E+00,
   /*   121  */  0.7804972E-06,
   /*   122  */  0.1167083E-05,
   /*   123  */  0.9855922E-06,
   /*   124  */  0.1377863E-06,
   /*   125  */  0.0000000E+00,
   /*   126  */  0.0000000E+00,
   /*   127  */  0.0000000E+00,
   /*   128  */  0.0000000E+00,
   /*   129  */  0.1054137E-04,
   /*   130  */  0.3000794E-04,
   /*   131  */  0.6023544E-04,
   /*   132  */  0.1054117E-03,
   /*   133  */  0.1708793E-03,
   /*   134  */  0.2687827E-03,
   /*   135  */  0.4137881E-03,
   /*   136  */  0.6199548E-03,
   /*   137  */  0.9161179E-03,
   /*   138  */  0.1342796E-02,
   /*   139  */  0.1943264E-02,
   /*   140  */  0.2764098E-02,
   /*   141  */  0.3827799E-02,
   /*   142  */  0.5169721E-02,
   /*   143  */  0.6821727E-02,
   /*   144  */  0.8854456E-02,
   /*   145  */  0.1132517E-01,
   /*   146  */  0.1430050E-01,
   /*   147  */  0.1780783E-01,
   /*   148  */  0.2188450E-01,
   /*   149  */  0.2658871E-01,
   /*   150  */  0.3198998E-01,
   /*   151  */  0.3815571E-01,
   /*   152  */  0.4513576E-01,
   /*   153  */  0.5298176E-01,
   /*   154  */  0.6175103E-01,
   /*   155  */  0.7152569E-01,
   /*   156  */  0.8236921E-01,
   /*   157  */  0.9432708E-01,
   /*   158  */  0.1074572E+00,
   /*   159  */  0.1218557E+00,
   /*   160  */  0.1376028E+00,
   /*   161  */  0.1547317E+00,
   /*   162  */  0.1733313E+00,
   /*   163  */  0.1934747E+00,
   /*   164  */  0.2152692E+00,
   /*   165  */  0.2387823E+00,
   /*   166  */  0.2641039E+00,
   /*   167  */  0.2913149E+00,
   /*   168  */  0.3204254E+00,
   /*   169  */  0.3516649E+00,
   /*   170  */  0.3850510E+00,
   /*   171  */  0.4207101E+00,
   /*   172  */  0.4586671E+00,
   /*   173  */  0.4990536E+00,
   /*   174  */  0.5420371E+00,
   /*   175  */  0.5876358E+00,
   /*   176  */  0.6359682E+00,
   /*   177  */  0.6837495E+00,
   /*   178  */  0.7288000E+00,
   /*   179  */  0.7716000E+00,
   /*   180  */  0.8113000E+00,
   /*   181  */  0.8474000E+00,
   /*   182  */  0.8797000E+00,
   /*   183  */  0.9079000E+00,
   /*   184  */  0.9319000E+00,
   /*   185  */  0.9518000E+00,
   /*   186  */  0.9676000E+00,
   /*   187  */  0.9797000E+00,
   /*   188  */  0.9883000E+00,
   /*   189  */  0.9940000E+00,
   /*   190  */  0.9976000E+00,
   /*   191  */  0.1000000E+01,
};
