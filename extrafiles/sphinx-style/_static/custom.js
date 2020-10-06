$(document).ready(function() {
  if ($('#fancy-particles-small').length) {
  	$.get('_static/particles.json')
    .done(function() { 
        particlesJS.load('fancy-particles-small', '_static/particles.json', null);
    }).fail(function() { 
        particlesJS.load('fancy-particles-small', '../_static/particles.json', null);
    })
  }
  // Change page content for cpp_api, since this is not possible with exhale
  var curUrl = $(location).attr("href");
  if(curUrl.indexOf('cpp_api') >= 0) {
    $('a[href$="#class-hierarchy"]').css('display','none');
    // Change sidebar content for cpp_api subpages, since this is not possible with exhale
    if(curUrl.indexOf('library_root') < 0) {
      $('#sidebar').html('<ul class="sidebar"><ul class="current"><li class="toctree-l1"><a class="reference internal" href="../python_api/modules.html">Python Documentation</a></li><li class="toctree-l1 current"><a class="current reference internal" href="library_root.html">C++ Documentation</a><ul><li class="toctree-l2"><a class="reference internal" href="library_root.html#file-hierarchy">File Hierarchy</a></li><li class="toctree-l2"><a class="reference internal" href="library_root.html#full-api">Full API</a><ul><li class="toctree-l3"><a class="reference internal" href="library_root.html#namespaces">Namespaces</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#classes-and-structs">Classes and Structs</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#enums">Enums</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#functions">Functions</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#variables">Variables</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#defines">Defines</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#typedefs">Typedefs</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#directories">Directories</a></li></ul></li></ul></li><li class="toctree-l1"><a class="reference internal" href="../notebooks.html">Jupyter Notebook</a></li><li class="toctree-l1"><a class="reference internal" href="../DevGuide.html">Developer Guide</a></li></ul></ul><form action="../search.html" method="get"><div class="form-group"><input type="text" name="q" class="form-control" placeholder="Search" /></div><input type="hidden" name="check_keywords" value="yes" /><input type="hidden" name="area" value="default" /></form>')
    }
  }
});