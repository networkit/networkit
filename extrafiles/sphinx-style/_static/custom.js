var docsJq = null;

function getDocsJq() {
  return docsJq || window.$jqTheme || window.jQuery || window.$;
}

function getDocsUrlRoot() {
  if (window.DOCUMENTATION_OPTIONS &&
      typeof window.DOCUMENTATION_OPTIONS.URL_ROOT === 'string' &&
      window.DOCUMENTATION_OPTIONS.URL_ROOT) {
    return window.DOCUMENTATION_OPTIONS.URL_ROOT;
  }

  // Sphinx 8+ exposes content root on the <html> tag.
  if (document.documentElement) {
    var contentRoot = document.documentElement.getAttribute('data-content_root');
    if (typeof contentRoot === 'string' && contentRoot) {
      return contentRoot;
    }
  }

  // Older Sphinx variants may provide it on the script tag.
  var optionsScript = document.getElementById('documentation_options');
  if (optionsScript) {
    var scriptRoot = optionsScript.getAttribute('data-url_root');
    if (typeof scriptRoot === 'string' && scriptRoot) {
      return scriptRoot;
    }
  }
  return './';
}

function normalizeDocsPath(urlPath) {
  var anchor = document.createElement('a');
  anchor.href = urlPath || '/';

  var path = anchor.pathname || urlPath || '/';
  if (path.charAt(0) !== '/') {
    path = '/' + path;
  }
  if (path.charAt(path.length - 1) !== '/') {
    path += '/';
  }
  return path;
}

function normalizeDocsHref(urlValue, baseHref, ensureTrailingSlash) {
  try {
    var url = new URL(urlValue || './', baseHref || window.location.href);
    if (ensureTrailingSlash !== false && url.pathname.charAt(url.pathname.length - 1) !== '/') {
      url.pathname += '/';
    }
    return url.href;
  } catch (error) {
    return null;
  }
}

function isRootRelativeUrl(urlValue) {
  return /^\/(?!\/)/.test(urlValue || '');
}

function pathLooksCompatible(pathValue) {
  var normalizedPath = normalizeDocsPath(pathValue);
  var currentPath = window.location.pathname || '/';
  if (currentPath.indexOf(normalizedPath) === 0) {
    return true;
  }

  if (/\/nightly\/$/.test(normalizedPath)) {
    var basePath = normalizedPath.replace(/nightly\/$/, '');
    if (currentPath.indexOf(basePath) === 0) {
      return true;
    }
  }

  if (/\/versions\.json$/.test(normalizedPath)) {
    var metadataBasePath = normalizedPath.replace(/versions\.json$/, '');
    if (currentPath.indexOf(metadataBasePath) === 0 ||
        currentPath.indexOf(metadataBasePath + 'nightly/') === 0) {
      return true;
    }
  }
  return false;
}

function resolveConfiguredDocsUrl(urlValue, fallbackUrl) {
  if (!urlValue) {
    return fallbackUrl;
  }

  if (isRootRelativeUrl(urlValue)) {
    if (window.location.protocol === 'file:' || !pathLooksCompatible(urlValue)) {
      return fallbackUrl;
    }
  }

  var normalizedHref = normalizeDocsHref(urlValue, window.location.href);
  return normalizedHref || fallbackUrl;
}

function resolveVersionsMetadataUrl(urlValue, stableUrl) {
  var fallbackUrl = normalizeDocsHref('versions.json', stableUrl, false);
  if (!urlValue) {
    return fallbackUrl;
  }

  if (isRootRelativeUrl(urlValue)) {
    if (window.location.protocol === 'file:' || !pathLooksCompatible(urlValue)) {
      return fallbackUrl;
    }
  }

  var normalizedHref = normalizeDocsHref(urlValue, window.location.href, false);
  return normalizedHref || fallbackUrl;
}

function detectCurrentDocsChannel(stableUrl, nightlyUrl, fallbackChannel) {
  var currentPath = window.location.pathname || '/';
  var nightlyPath = normalizeDocsPath(nightlyUrl);
  var stablePath = normalizeDocsPath(stableUrl);

  if (currentPath.indexOf(nightlyPath) === 0) {
    return 'nightly';
  }
  if (currentPath.indexOf(stablePath) === 0) {
    return 'stable';
  }
  return fallbackChannel || 'stable';
}

function getCurrentRelativePath(rootUrl) {
  var currentPath = window.location.pathname || '/';
  var rootPath = normalizeDocsPath(rootUrl);
  if (currentPath.indexOf(rootPath) === 0) {
    return currentPath.substring(rootPath.length);
  }

  var lastSlash = currentPath.lastIndexOf('/');
  if (lastSlash >= 0) {
    return currentPath.substring(lastSlash + 1);
  }
  return '';
}

function populateVersionSelect($select, config, currentChannel) {
  var jq = getDocsJq();
  if (!jq) {
    return;
  }
  var entries = ['stable', 'nightly'];
  $select.empty();

  for (var i = 0; i < entries.length; i += 1) {
    var channel = entries[i];
    if (!config[channel] || !config[channel].label) {
      continue;
    }
    var option = jq('<option></option>');
    option.attr('value', channel);
    option.text(config[channel].label);
    $select.append(option);
  }

  if ($select.find('option[value="' + currentChannel + '"]').length) {
    $select.val(currentChannel);
  } else if ($select.find('option[value="stable"]').length) {
    $select.val('stable');
  }
}

function setupDocsVersionSwitcher() {
  var jq = getDocsJq();
  if (!jq) {
    return;
  }

  var $select = jq('#docs-version-select');
  if (!$select.length) {
    return;
  }

  var docsRootUrl = normalizeDocsHref(getDocsUrlRoot(), window.location.href) ||
    normalizeDocsHref('./', window.location.href);
  var currentChannel = /\/nightly\/$/.test(normalizeDocsPath(docsRootUrl)) ? 'nightly' : 'stable';
  if ($select.attr('data-docs-channel')) {
    currentChannel = $select.attr('data-docs-channel');
  }

  var stableUrl = /\/nightly\/$/.test(normalizeDocsPath(docsRootUrl)) ?
    normalizeDocsHref('../', docsRootUrl) :
    docsRootUrl;
  var nightlyUrl = /\/nightly\/$/.test(normalizeDocsPath(docsRootUrl)) ?
    docsRootUrl :
    normalizeDocsHref('nightly/', docsRootUrl);

  stableUrl = resolveConfiguredDocsUrl($select.attr('data-docs-stable-url'), stableUrl);
  nightlyUrl = resolveConfiguredDocsUrl($select.attr('data-docs-nightly-url'), nightlyUrl);

  var config = {
    stable: {
      label: $select.attr('data-docs-stable-label') || 'Stable',
      url: stableUrl
    },
    nightly: {
      label: $select.attr('data-docs-nightly-label') || 'Nightly',
      url: nightlyUrl
    }
  };

  var currentLabel = $select.attr('data-docs-version-label');
  if (currentLabel) {
    if (currentChannel === 'stable' && !$select.attr('data-docs-stable-label')) {
      config.stable.label = currentLabel;
    }
    if (currentChannel === 'nightly' && !$select.attr('data-docs-nightly-label')) {
      config.nightly.label = currentLabel;
    }
  }

  currentChannel = detectCurrentDocsChannel(config.stable.url, config.nightly.url, currentChannel);
  populateVersionSelect($select, config, currentChannel);

  $select.off('change.docsVersion').on('change.docsVersion', function() {
    var targetChannel = $select.val();
    var targetBaseUrl = targetChannel === 'nightly' ? config.nightly.url : config.stable.url;
    var relativePath = getCurrentRelativePath(docsRootUrl);
    var targetUrl = new URL(relativePath, targetBaseUrl);
    targetUrl.search = window.location.search;
    targetUrl.hash = window.location.hash;
    window.location.assign(targetUrl.href);
  });

  var versionsUrl = resolveVersionsMetadataUrl(
    $select.attr('data-docs-versions-url'),
    config.stable.url
  );
  jq.getJSON(versionsUrl).done(function(metadata) {
    if (!metadata || !metadata.channels) {
      return;
    }

    if (metadata.channels.stable && metadata.channels.stable.url) {
      config.stable.url = resolveConfiguredDocsUrl(metadata.channels.stable.url, config.stable.url);
    }
    if (metadata.channels.nightly && metadata.channels.nightly.url) {
      config.nightly.url = resolveConfiguredDocsUrl(metadata.channels.nightly.url, config.nightly.url);
    }

    if (metadata.channels.stable && metadata.channels.stable.label) {
      config.stable.label = metadata.channels.stable.label;
    }
    if (metadata.channels.nightly && metadata.channels.nightly.label) {
      config.nightly.label = metadata.channels.nightly.label;
    }

    currentChannel = detectCurrentDocsChannel(config.stable.url, config.nightly.url, currentChannel);
    populateVersionSelect($select, config, currentChannel);
  });
}

document.addEventListener('DOMContentLoaded', function() {
  var jq = getDocsJq();
  if (!jq) {
    return;
  }
  docsJq = jq;

  var docsVersionSwitcherHtml = '';
  if (jq('#docs-version-switcher').length) {
    docsVersionSwitcherHtml = jq('#docs-version-switcher').prop('outerHTML');
  }

  if (jq('#fancy-particles-small').length) {
  	jq.get('_static/particles.json')
    .done(function() { 
        particlesJS.load('fancy-particles-small', '_static/particles.json', null);
    }).fail(function() { 
        particlesJS.load('fancy-particles-small', '../_static/particles.json', null);
    })
  }
  // Change page content for cpp_api, since this is not possible with exhale
  var curUrl = jq(location).attr("href");
  if(curUrl.indexOf('cpp_api') >= 0) {
    jq('a[href$="#class-hierarchy"]').css('display','none');
    // Change sidebar content for cpp_api subpages, since this is not possible with exhale
    if(curUrl.indexOf('library_root') < 0) {
      jq('#sidebar').html('<ul class="sidebar"><ul class="current"><li class="toctree-l1"><a class="reference internal" href="../python_api/modules.html">Python Documentation</a></li><li class="toctree-l1 current"><a class="current reference internal" href="library_root.html">C++ Documentation</a><ul><li class="toctree-l2"><a class="reference internal" href="library_root.html#file-hierarchy">File Hierarchy</a></li><li class="toctree-l2"><a class="reference internal" href="library_root.html#full-api">Full API</a><ul><li class="toctree-l3"><a class="reference internal" href="library_root.html#namespaces">Namespaces</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#classes-and-structs">Classes and Structs</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#enums">Enums</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#functions">Functions</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#variables">Variables</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#defines">Defines</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#typedefs">Typedefs</a></li><li class="toctree-l3"><a class="reference internal" href="library_root.html#directories">Directories</a></li></ul></li></ul></li><li class="toctree-l1"><a class="reference internal" href="../notebooks.html">Jupyter Notebook</a></li><li class="toctree-l1"><a class="reference internal" href="../DevGuide.html">Developer Guide</a></li></ul></ul><form action="../search.html" method="get"><div class="form-group"><input type="text" name="q" class="form-control" placeholder="Search" /></div><input type="hidden" name="check_keywords" value="yes" /><input type="hidden" name="area" value="default" /></form>')
      if (docsVersionSwitcherHtml) {
        jq('#sidebar').prepend(docsVersionSwitcherHtml);
      }
    }
  }

  setupDocsVersionSwitcher();
});
