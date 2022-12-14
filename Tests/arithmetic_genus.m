printf "Testing arithmetic genus... D=";
arithmetic_genus_table :=
[ <5, 1>,
<8, 1>,
<12, 1>,
<13, 1>,
<17, 1>,
<21, 1>,
<24, 1>,
<28, 1>,
<29, 2>,
<33, 1>,
<37, 2>,
<40, 2>,
<41, 2>,
<44, 2>,
<53, 3>,
<56, 2>,
<57, 2>,
<60, 1>,
<61, 3>,
<65, 3>,
<69, 2>,
<73, 3>,
<76, 3>,
<77, 3>,
<85, 4>,
<88, 3>,
<89, 4>,
<92, 3>,
<93, 3>,
<97, 4>,
<101, 5>,
<104, 5>,
<105, 2>,
<109, 5>,
<113, 5>,
<120, 3>,
<124, 4>,
<129, 6>,
<133, 5>,
<136, 2>,
<137, 6>,
<140, 4>,
<141, 4>,
<145, 7>,
<149, 7>,
<152, 6>,
<156, 5>,
<157, 7>,
<161, 5>,
<165, 4>,
<168, 4>,
<172, 8>,
<173, 8>,
<177, 7>,
<181, 8>,
<184, 6>,
<185, 9>,
<188, 6>,
<193, 10>,
<197, 9>,
<201, 10>,
<204, 8>,
<205, 6>,
<209, 10>,
<213, 6>,
<217, 9>,
<220, 7>,
<221, 6>,
<229, 10>,
<232, 11>,
<233, 12>,
<236, 10>,
<237, 9>,
<241, 14>,
<248, 9>,
<249, 13>,
<253, 9>,
<257, 13>,
<264, 11>,
<265, 15>,
<268, 13>,
<269, 12>,
<273, 10>,
<277, 14>,
<280, 10>,
<281, 16>,
<284, 11>,
<285, 8>,
<293, 13>,
<296, 16>,
<301, 13>,
<305, 13>,
<309, 13>,
<312, 9>,
<313, 19>,
<316, 8>,
<317, 14>,
<321, 11>,
<328, 17>,
<329, 15>,
<332, 15>,
<337, 22>,
<341, 13>,
<344, 17>,
<345, 14>,
<348, 13>,
<349, 17>,
<353, 19>,
<357, 10>,
<364, 18>,
<365, 18>,
<373, 20>,
<376, 17>,
<377, 17>,
<380, 15>,
<381, 17>,
<385, 20>,
<389, 19>,
<393, 23>,
<397, 19>,
<401, 25>,
<408, 18>,
<409, 29>,
<412, 21>,
<413, 17>,
<417, 26>,
<421, 22>,
<424, 26>,
<428, 22>,
<429, 14>,
<433, 30>,
<437, 17>,
<440, 17>,
<444, 21>,
<445, 24>,
<449, 29>,
<453, 19>,
<456, 23>,
<457, 33>,
<460, 20>,
<461, 22>,
<465, 24>,
<469, 17>,
<472, 25>,
<473, 20>,
<476, 18>,
<481, 37>,
<485, 24>,
<488, 25>,
<489, 34>,
<492, 24>,
<493, 26>,
<497, 25> ];


for _ in [1..10] do
  D, chi := Explode(Random(arithmetic_genus_table));
  printf "%o ", D;
  assert chi eq ArithmeticGenus(QuadraticField(D));
end for;