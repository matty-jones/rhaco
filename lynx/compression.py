import zlib
import base64

a = '''     ,  ,\n    (\ "\n    ,--;.)._\n   ).,-._ . ""-,_\n  /.'".- " 8 o . ";_                             \n . 8.o .""-.---...,,--------.._   _"";\n   """  ")) 8 . . 8 . 8   8  8  8  8. 8 8 ._""._;\n         ";. .8 .8  .8  8  8  8  8 . 8. 8 .".""\n            ;.. 8 ; .  8. 8  8  8 . } 8 . 8 :\n             ;.. 8 ; 8. 8  8  8  8 (  . 8 . :\n               ;. 8 \ .   .......;;;  8 . 8 :\n                ;o  ;"\\\\```````( o(  8   .;\n                : o:  ;           :. : . 8 (\n                :o ; ;             "; ";. o :\n                ; o; ;               "; ;";..\ \n                ;.; .:                )./  ;. ;\n               _).< .;              _;./  _;./\n             ;"__/--"             ((__7  ((_J'''


encoded = base64.encodestring(zlib.compress(a))
print(encoded)
decoded = zlib.decompress(base64.decodestring(encoded))
#print(decoded)

# pyLZMA
