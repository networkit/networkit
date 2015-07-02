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


function showOverlay(source)
{
	var images = source.getAttribute("data-image");
	var element = document.getElementById("NetworKit_Overlay_Image");
	element.setAttribute("data-image", images);
	element.style.backgroundImage = "url(" + images +")";
	showElement("NetworKit_Overlay", true);
}


function showElement(id, show)
{
	var element = document.getElementById(id);
	element.style.display = (show) ? "block" : "none";
}