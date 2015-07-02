function embedPlots(id)
{
	var elements = document.getElementById(id).getElementsByClassName("Plot");
	var i;
	for (i=0; i<elements.length; i++) {
		var content =
			"<div class=\\"Image\\" style=\\"background-image:url(" + elements[i].getAttribute("data-image") + ")\\">" +
			"  <div class=\\"Title\\">" + elements[i].getAttribute("title") + "</div>" +
			"</div>";
		elements[i].innerHTML = content;
		elements[i].onclick = function (e) {
			showOverlay((e.target) ? e.target : e.srcElement);
		}
	}
}


function showElement(id, show)
{
	var element = document.getElementById(id);
	element.style.display = (show) ? "block" : "none";
}


function showOverlay(source)
{
	var data = source.getAttribute("data-image");
	var image = document.getElementById("NetworKit_Overlay_Image");
	image.setAttribute("data-image", data);
	image.style.backgroundImage = "url(" + data +")";
	var link = document.getElementById("NetworKit_Overlay_Toolbar_Save_Link");
	link.href = data;
	link.download = "image.svg";
	showElement("NetworKit_Overlay", true);
}
