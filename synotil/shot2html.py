import synotil.dio as sdio
# import PIL
# import cairosvg
# import pymisca.util as pyutil
import os
def make_thumbnail(path, size):
#     newpath = os.path.dirname(path)
#     image_name = os.path.basename(os.path.splitext(path)[0])
#     newpath = newpath + '/' + image_name + '.thumbnail.jpg'


    try:
        sp = path.rsplit('.',1)
        if len(sp)==1:
            ext = ''
            base = sp[0]
        else:
            base,ext = sp
        if ext.lower()=='svg':
            base = path.rsplit('.',1)[0]
            with open(path,'r') as f:
                svg = f.read()
            cairosvg.svg2png(bytestring=svg,write_to=base + '.png')
            path = base + '.png'
            ext = 'png'

        newpath =  base + '__thumbnail.' + ext 
        img = PIL.Image.open(path)
        img.thumbnail(size)
        img.save(newpath, None)

    except Exception as e:
        newpath = path
        print('[WARN]image error for PATH %s:\n%s'%(newpath,e))
    return newpath


def shot2html(indexFile,localPath = True):
# indexFile = "PROG=chipShots_bedFile=64-SD-ZT20-ELF3OX-RERUN_S8_peaks.narrowPeak.FC-GT-3dot0/index.tsv"
    dfc = sdio.extract_peak(indexFile)
    dfc['img'] = dfc[dfc.columns[-1]]
#     THUMBNAIL_SIZE = 256, 256
    # HTML_REPORT_NAME='test.html'
    ofname = HTML_REPORT_NAME='%s.html'%indexFile 
    TITLE = os.path.basename(indexFile)
    # !pip install --user cairosvg==1.0.22
    # cairosvg.svg2png?
    # from cairosvg import svg2png

    # generate an html page

    fhtml = open(HTML_REPORT_NAME, 'w')
#     TITLE = 'test'
    fhtml.write('''
    <html>
    <head>
    <style> 
        body {{font-family:\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica neue\"}}
    </style>
    <title>IGV screenshots for {TITLE}</title>
    '''.format(**locals()))

    js = """
    <script>
    //var frame = window.frames["postFrame"];
    //var frame = document.getElementsByName("postFrame")[0];

    //dir = document.URL.replace(/^.*[\\\/]/, '')
    //var frame = window.open("about:blank","postFrame", 'height=440,width=260,scrollbars=yes');
    console.log(document.URL);


    // Pass the checkbox name to the function
    function checkbox() {
       // var frame = window.frames["postFrame"];
    //    var doc = frame.document;
    //  var doc = frame.contentDocument || frame.contentWindow.document

      //opened.open();
      var checkboxes = document.getElementsByName("line");
      var checkboxesChecked = [];

      // loop over them all
      html_string = '';
      for (var i=0; i<checkboxes.length; i++) {
         // And stick the checked ones onto an array...
         if (checkboxes[i].checked) {
            checkboxesChecked.push(checkboxes[i]);
            html_string += checkboxes[i].value + '<br>';
    //        doc.writeln(checkboxes[i].value);
    //        doc.write("<br>");
         }
      }
      src = "data:text/html;charset=utf-8," + escape(html_string);
       var frame = window.open(src,"postFrame", 'height=440,width=260,scrollbars=yes');

      //opened.close()
      // Return the array if it is non-empty, or null
      return checkboxesChecked.length > 0 ? checkboxesChecked : null;
    }
    </script>
    </head>

    """

    fhtml.write(js)



    fhtml.write('<body>\n')

    # fc = []
    # for t in image_lst:
    #     narrowpeak_line = t[0]
    #     lst = narrowpeak_line.split()
    #     fc.append(float(lst[4]))

    # max_fc = max(fc)
    # min_fc = min(fc)
    # lmd = 0.8
    # show_color_val = (1.0 - lmd) * min_fc +  1.0 * lmd * max_fc


    colors = ['#ff9999','#ff8080','#ff6666','#ff4d4d','#ff3333','#ff1a1a','#ff0000']
    sz_color = len(colors)

    fhtml.write('<p><b>%s</b></p>' % ('&nbsp; &nbsp;'.join(['Peak_position','Peak_name', 'Integer_score_for_display', 'Dot','fold-change',' -log10pvalue',' -log10qvalue','Relative_summit_position_to_peak_start']))) # head

    fhtml.write('<iframe name="postFrame"></iframe>')
    fhtml.write('<form>')
    fhtml.write('<input type=\"submit\" value=\"Save\" size=\"35\" onClick=\"return checkbox();\">\n')

    for i in range(len(dfc)):
        row = dfc.iloc[i]
    #     if np.isnan(row.img):

        if isinstance(row.img,float):
            imgFile = './none'
        else:
            imgFile = row.img
            
        if localPath:
            imgFile = './%s'% os.path.basename(imgFile)
        else:
            imgFile = '/notebooks%s' % imgFile            
            
        bedLine = ('\t').join(map(str,row[:-1]))
    #     bedLine = '\t'.join(row)
    #     line = line_lst[i]
        ind = i
        fhtml.write('<input type=\"checkbox\" id=\"%s\" name=\"line\" value=\"%s\" >' % (ind, 
                                                                                         bedLine.replace('\t','\t')
                                                                                        ))
        c = colors[0]
        fhtml.write('<font color=%s>%s</font><br/>\n' % (c,bedLine))
        fhtml.write('<font color=%s>%s</font><br/>\n' % (c, bedLine.replace('\t','&emsp; ')))
    #     timgFile = make_thumbnail(imgFile, THUMBNAIL_SIZE)
    #     if JUPYTER:

    
    #     timgFile = '/notebooks' + timgFile
    #     if imgFile.endswith('.svg'):
    #         imgFile = imgFile.replace('.svg','.png')

        fhtml.write('<a href=\"%s\">\n<img width=800 src=\"%s\"></a><br/><br/>\n' % (imgFile, imgFile))

    fhtml.write('</form>')    
    fhtml.write('</body>\n')
    fhtml.write('</html>\n')
    fhtml.close()    

    # os.remove(PYTHON_SCRIPT_NAME)

    print('[Done.]')
    return ofname
