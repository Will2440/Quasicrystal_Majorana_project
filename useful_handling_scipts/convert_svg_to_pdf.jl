using Base.Filesystem
using Cairo
using Rsvg

function convert_svg_to_pdf(svg_path, pdf_path)
    try
        handle = Rsvg.handle_new_from_file(svg_path)
        dims = Rsvg.handle_get_dimensions(handle)
        width = max(dims.width, 1)
        height = max(dims.height, 1)
        surface = Cairo.CairoPDFSurface(pdf_path, width, height)
        ctx = Cairo.CairoContext(surface)
        Rsvg.handle_render_cairo(handle, ctx.ptr)
        Cairo.finish(surface)
        Cairo.destroy(surface)
        println("Converted $svg_path to $pdf_path")
    catch e
        println("Error converting $svg_path: $e")
    end
end

function batch_convert_svgs_to_pdfs(directory)
    for filename in readdir(directory)
        if endswith(filename, ".svg")
            svg_path = joinpath(directory, filename)
            pdf_dir = joinpath(directory, "pdf_versions")
            if !isdir(pdf_dir)
                mkdir(pdf_dir)
            end
            pdf_path = joinpath(pdf_dir, replace(filename, ".svg" => ".pdf"))
            convert_svg_to_pdf(svg_path, pdf_path)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    current_dir = pwd()
    println("Converting SVGs in: $current_dir")
    batch_convert_svgs_to_pdfs(current_dir)
end