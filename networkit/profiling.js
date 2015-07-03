function pageEmbed(id)
{
	var elements = document.getElementById(id).getElementsByClassName("Plot");
	var i, j;
	for (i=0; i<elements.length; i++) {
		elements[i].id = id + "_Plot_" + i;
		var data = elements[i].getAttribute("data-image").split("|");
		elements[i].removeAttribute("data-image");
		var content = 
			"<div class=\\"Image\\" id=\\"" + elements[i].id + "_Image\\">" +
			"  <div class=\\"Title\\">" + elements[i].getAttribute("title") + "</div>" +
			"</div>";
		elements[i].innerHTML = content;
		elements[i].setAttribute("data-image-index", 0);
		elements[i].setAttribute("data-image-length", data.length);
		for (j=0; j<data.length; j++) {
			elements[i].setAttribute("data-image-" + j, data[j]);
		}
		plotUpdate(elements[i]);
		elements[i].onclick = function (e) {
			overlayShow((e.target) ? e.target : e.srcElement);
		}
	}
}


function plotUpdate(source)
{
	var index = source.getAttribute("data-image-index");
	var data = source.getAttribute("data-image-" + index);
	var image = document.getElementById(source.id + "_Image");
	image.style.backgroundImage = "url(" + data + ")";
}


function showElement(id, show)
{
	var element = document.getElementById(id);
	element.style.display = (show) ? "block" : "none";
}


function overlayShow(source)
{
	overlayUpdate(source);
	showElement("NetworKit_Overlay", true);
}


function overlayUpdate(source)
{
	var index = source.getAttribute("data-image-index");
	var data = source.getAttribute("data-image-" + index);
	var image = document.getElementById("NetworKit_Overlay_Image");
	image.setAttribute("data-id", source.id);
	image.style.backgroundImage = "url(" + data + ")";
	var link = document.getElementById("NetworKit_Overlay_Toolbar_Save_Link");
	link.href = data;
	link.download = "image.svg";
}


function overlayImageShift(delta)
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
	overlayUpdate(source);
}