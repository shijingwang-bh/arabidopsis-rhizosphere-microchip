"""Assemble current Figure 1, Figure 2, S1 and S4 from preserved source assets."""
from pathlib import Path
from io import BytesIO
import shutil
import sys
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas

pkg=Path(sys.argv[1]).resolve();run=Path(sys.argv[2]).resolve()
fig=run/'figures';parts=run/'parts';fig.mkdir(parents=True,exist_ok=True)
mm=72/25.4
shutil.copy2(pkg/'assets/platform/design_and_hardware.pdf',fig/'Figure_1.pdf')
shutil.copy2(pkg/'assets/microscopy/longitudinal_plate.pdf',fig/'Figure_S4.pdf')

# Retain the original R vector lateral-root and eight-gene panels.
writer=PdfWriter();page=writer.add_blank_page(width=186*mm,height=78*mm)
for name,x,width in [('lateral_roots.pdf',0,64),('deg_effects.pdf',69,116)]:
 panel=PdfReader(parts/name).pages[0]
 page.merge_transformed_page(panel,Transformation().scale(width*mm/float(panel.mediabox.width)).translate(x*mm,0))
buf=BytesIO();c=canvas.Canvas(buf,pagesize=(186*mm,78*mm))
for letter,title,x in [('a','Lateral-root number',1),('b','Leaf differential expression',70)]:
 c.setFont('Helvetica-Bold',12);c.drawString(x*mm,74*mm,letter)
 c.setFont('Helvetica-Bold',10);c.drawString((x+6)*mm,74*mm,title)
c.save();buf.seek(0);page.merge_page(PdfReader(buf).pages[0])
bottom=PdfReader(fig/'Figure_2_bottom.pdf').pages[0]
top_h=float(page.mediabox.height);bottom_h=float(bottom.mediabox.height);gap=13
combined_width=float(bottom.mediabox.width)
combined=PdfWriter();target=combined.add_blank_page(width=combined_width,height=top_h+bottom_h+gap)
target.merge_page(bottom)
target.merge_transformed_page(page,Transformation().translate((combined_width-float(page.mediabox.width))/2,bottom_h+gap))
with (fig/'Figure_2.pdf').open('wb') as f:combined.write(f)

# The original a-c image plate is preserved; d is generated from all four
# pre-inoculation BH tests, then only the three cultivation traits are plotted.
asset=pkg/'assets/platform/S1_panels_abc.png'
plot=PdfReader(parts/'cultivation.pdf').pages[0]
width=190*mm
with Image.open(asset) as im:plate_h=width*im.height/im.width
plot_w=490.;scale=plot_w/float(plot.mediabox.width);plot_h=float(plot.mediabox.height)*scale
height=plate_h+12+18+plot_h+5
buf=BytesIO();c=canvas.Canvas(buf,pagesize=(width,height))
c.drawImage(str(asset),0,height-plate_h,width=width,height=plate_h)
c.setFont('Helvetica-Bold',13);c.drawString(6,plot_h+16,'d')
c.setFont('Helvetica-Bold',10.5);c.drawString(27,plot_h+16,'Cultivation performance')
c.save();buf.seek(0)
writer=PdfWriter();page=writer.add_page(PdfReader(buf).pages[0])
page.merge_transformed_page(plot,Transformation().scale(scale).translate((width-plot_w)/2,5))
with (fig/'Figure_S1.pdf').open('wb') as f:writer.write(f)
print('Assembled current Figures 1, 2, S1 and S4')
