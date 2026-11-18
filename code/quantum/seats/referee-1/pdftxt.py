import sys
import pymupdf
d = pymupdf.open(sys.argv[1])
open(sys.argv[2], 'w', encoding='utf-8').write('\n'.join(p.get_text() for p in d))
