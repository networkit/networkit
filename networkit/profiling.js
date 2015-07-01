function embedPlot(id,
                   width)
{
	var elements = document.getElementById(id).getElementsByClassName("Plot");
	var i;
	for (i=0; i<elements.length; i++) {
		elements[i].innerHTML = "<img width=\\"" + width + "\\" src=\\"" + elements[i].getAttribute("data-image") + "\\" />";
		elements[i].setAttribute("data-image", "");
	}
}