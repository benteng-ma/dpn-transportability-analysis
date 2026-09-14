from pathlib import Path
from PIL import Image,ImageOps,ImageDraw
B=Path(__file__).resolve().parents[1]
pages=sorted((B/'rendered').glob('page-*.png'))
for start in range(0,len(pages),4):
 canvas=Image.new('RGB',(1200,1700),'white');d=ImageDraw.Draw(canvas)
 for j,p in enumerate(pages[start:start+4]):
  im=Image.open(p).convert('RGB');im.thumbnail((590,800));x=(j%2)*600;y=(j//2)*850;canvas.paste(im,(x,y+25));d.text((x+10,y+3),p.stem,fill='black')
 canvas.save(B/'rendered'/f'contact-{start//4+1}.png')
