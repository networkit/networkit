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
  }
});