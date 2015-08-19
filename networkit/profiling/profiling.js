/*
	file: profiling.js
	author: Mark Erb
		
	This file contains all JavaScript scripts for the profiling layout 
	
	The data will be automatically embedded as singelton into the original
	IPython Layout 
	
		html > head > script id="NetworKit_script"
		
	In addition a function for hiding the Overlay is defined in 
	"profiling.py"
	
	To prevent conflicts with the IPython Notebook, prove that every function
	begins with	"NetworKit_".
*/

/*
	Markup the data from "profiling.profile.html"
	
	Arguments:
		id: ID of the "NetworKit_Page" to markup
 */
function NetworKit_pageEmbed(id)
{
	var i, j;
	var elements;
	
	elements = document.getElementById(id).getElementsByClassName("Plot");
	for (i=0; i<elements.length; i++) {
		elements[i].id = id + "_Plot_" + i;
		var data = elements[i].getAttribute("data-image").split("|");
		elements[i].removeAttribute("data-image");
		var content = 
			"<div class=\\"Image\\" id=\\"" + elements[i].id + "_Image\\" />";
		elements[i].innerHTML = content;
		elements[i].setAttribute("data-image-index", 0);
		elements[i].setAttribute("data-image-length", data.length);
		for (j=0; j<data.length; j++) {
			elements[i].setAttribute("data-image-" + j, data[j]);
		}
		NetworKit_plotUpdate(elements[i]);
		elements[i].onclick = function (e) {
			NetworKit_overlayShow((e.target) ? e.target : e.srcElement);
		}
	}
	
	elements = document.getElementById(id).getElementsByClassName("HeatCell");
	for (i=0; i<elements.length; i++) {
		var data;
		
		data = Math.abs(parseFloat(elements[i].getAttribute("data-heat")));
		elements[i].style.backgroundColor = (data <= 1) ? "hsl(" + (240 + 120 * data) + ", 60%, 70%)" : "#00FF00";
		
		data = elements[i].getAttribute("data-image").split("|");
		elements[i].removeAttribute("data-image");
		elements[i].setAttribute("data-image-index", 0);
		elements[i].setAttribute("data-image-length", data.length);
		for (j=0; j<data.length; j++) {
			elements[i].setAttribute("data-image-" + j, data[j]);
		}
		elements[i].onclick = function (e) {
			NetworKit_overlayShow((e.target) ? e.target : e.srcElement);
		}
	}
	
	elements = document.getElementById(id).getElementsByClassName("Details");
	for (i=0; i<elements.length; i++) {
		elements[i].setAttribute("data-title", "-");
		NetworKit_toggleDetails(elements[i]);
		elements[i].onclick = function (e) {
			NetworKit_toggleDetails((e.target) ? e.target : e.srcElement);
		}
	}
	
	elements = document.getElementById(id).getElementsByClassName("MathValue");
	for (i=elements.length-1; i>=0; i--) {
		value = elements[i].innerHTML.trim();
		if (value === "nan") {
			elements[i].parentNode.innerHTML = ""
		}
	}
	
	elements = document.getElementById(id).getElementsByClassName("SubCategory");
	for (i=elements.length-1; i>=0; i--) {
		value = elements[i].innerHTML.trim();
		if (value === "") {
			elements[i].parentNode.removeChild(elements[i])
		}
	}
	
	elements = document.getElementById(id).getElementsByClassName("Category");
	for (i=elements.length-1; i>=0; i--) {
		value = elements[i].innerHTML.trim();
		if (value === "") {
			elements[i].parentNode.removeChild(elements[i])
		}
	}
	
	var isFirefox = false;
	try {
		isFirefox = typeof InstallTrigger !== "undefined";
	}
	catch (e) {}
	if (!isFirefox) {
		alert("Currently the function\'s output is only fully supported by Firefox.");
	}
}


/*
	Update Thumbnails.
	
	Arguments:
		source: Element (class=NetworKit_Plot) to update
*/
function NetworKit_plotUpdate(source)
{
	var index = source.getAttribute("data-image-index");
	var data = source.getAttribute("data-image-" + index);
	var image = document.getElementById(source.id + "_Image");
	image.style.backgroundImage = "url(" + data + ")";
}


/*
	"Hide/show" an element.
	
	Arguments:
		id: ID of the element to hide/show
		show: true  -> show
			  false -> hide
*/
function NetworKit_showElement(id, show)
{
	var element = document.getElementById(id);
	element.style.display = (show) ? "block" : "none";
}


/*
	Show Overlay.
	
	Arguments:
		source: element to get data for the overlay
*/
function NetworKit_overlayShow(source)
{
	NetworKit_overlayUpdate(source);
	NetworKit_showElement("NetworKit_Overlay", true);
}


/*
	Update data of overlay.

	Arguments:
		source: element to get data for the overlay
*/
function NetworKit_overlayUpdate(source)
{
	document.getElementById("NetworKit_Overlay_Title").innerHTML = source.title;
	var index = source.getAttribute("data-image-index");
	var data = source.getAttribute("data-image-" + index);
	var image = document.getElementById("NetworKit_Overlay_Image");
	image.setAttribute("data-id", source.id);
	image.style.backgroundImage = "url(" + data + ")";
	var link = document.getElementById("NetworKit_Overlay_Toolbar_Bottom_Save");
	link.href = data;
	link.download = source.title + ".svg";
}


/*
	Shift through the possible images of the overlay.
	
	Arguments:
		delta: >0 next delta-th picture
		       <0 previous delta-th picture
*/
function NetworKit_overlayImageShift(delta)
{
	var image = document.getElementById("NetworKit_Overlay_Image");
	var source = document.getElementById(image.getAttribute("data-id"));
	var index = parseInt(source.getAttribute("data-image-index"));
	var length = parseInt(source.getAttribute("data-image-length"));
	var index = (index+delta) % length;
	if (index < 0) {
		index = length + index;
	}
	source.setAttribute("data-image-index", index);
	NetworKit_overlayUpdate(source);
}


/*
	Toggle Measure Details
	
	Arguments:
		source: element to get data for the overlay
*/
function NetworKit_toggleDetails(source)
{
	var childs = source.children;
	var show = false;
	if (source.getAttribute("data-title") == "-") {
		source.setAttribute("data-title", "+");
		show = false;
	}
	else {
		source.setAttribute("data-title", "-");
		show = true;
	}
	for (i=0; i<childs.length; i++) {
		if (show) {
			childs[i].style.display = "block";
		}
		else {
			childs[i].style.display = "none";
		}
	}
}